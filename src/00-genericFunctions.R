
# Functions that might be used in another analysis

# get shared coordinates (only for single base)
# for example get all CpGs that are present in all samples
# @param coordDTList a list with each item being a data.frame for a separate 
# sample. The data frame must have a column chrBase that has chromosome
# and start information in a single string

getSharedCoord <- function(coordDTList) {
    
    # I first get vector of all cpgs (found in at least one sample)
    # because if they are not in the first sample, they are not in all samples 
    # so just take cpgs in first sample
    sharedCVec = coordDTList[[1]]$chrBase
    
    # progressively eliminate cpgs that are not in other samples
    for (i in 2:length(coordDTList)) {
        
        # only keeps shared cytosines each time
        sharedCVec = sharedCVec[sharedCVec %in% coordDTList[[i]]$chrBase]
        
    }
    
    # cytosine coordinate was not always in order of increasing coordinate
    sharedCVec = sort(sharedCVec, decreasing = FALSE)
    
    return(sharedCVec)
}


consensusMethyl = function(x, y) {
    # selecting consensus cytosines then ordering in ascending order
    # according to chrBase column
    data.table::setorder(x[x$chrBase %in% y, ], chrBase)
}

# filter based on input set of coordinates
# by reference?
filterCoord <- function(coordDTList, sharedCVec) {
    sharedCDTList = lapply(coordDTList, consensusMethyl, y=sharedCVec)
    lapply(sharedCDTList, function(x) x[, chrBase := NULL])
    return(sharedCDTList)
}


# validate that objects are matching dims and names
validateObjects = function(obj1, obj2, valDim=TRUE, valNames=FALSE) {
    if (valDim) {   
        if (is.null(dim(obj1))) {
            dim1 = length(obj1)
        } else {
            dim1 = dim(obj1)
        }
        
        if (is.null(dim(obj2))) {
            dim2 = length(obj2)
        } else {
            dim2 = dim(obj2)
        }
        
        # if either dim is 1, make sure larger dim matches other object
        # will break if both dimensions are 1
        if (any(dim1 %in% 1)) {
            dim1 = dim1[!(dim1 %in% 1)]
        }
        if (any(dim2 %in% 1)) {
            dim2 = dim2[!(dim2 %in% 1)]
        }
        
        # now check whether dims match
        if (!(any(dim2 %in% dim1) && any(dim1 %in% dim2))) {
            stop("Object dimensions do not match.")
        }
    }
    
}


# general pipeline/workflow for MIRA so I don't have to keep copying code
# required inputs are methFiles-RRBS files with methylation calls and
# region sets, MIRA package must already be loaded
# NOTE: Resize regions beforehand
# @param BSDTList List of data.tables compatible with aggregateMethyl.
# @param regionGRList GRangesList of region sets. Should be named but 
# will be given automatic names.
# IMPORTANT: resize beforehand so all regions are same size and appropriate length
# @param binNum Number of bins for aggregateMethyl()
# @param minBaseCovPerBin See ?aggregateMethyl
# @param annoDT 
# @param methylNormalizationVal Binned methylation data will be normalized by
# subtracting the minimum value from each MIRA profile (so that the lowest
# value for all profiles will be the same and the shape may be more easily 
# compared). Then methylNormalizationVal is added to all every bin so that the lowest bin
# is no longer zero but is this param instead (prevents division by zero).

MIRAPipeline = function(BSDTList, regionGRList, binNum=21, 
                        minBaseCovPerBin, annoDT, methylNormalizationVal=0.05) {
    # making sure arguments with no default are given
    if (any(c(missing(BSDTList), missing(regionGRList), 
              missing(minBaseCovPerBin), missing(annoDT)))) {
        stop("Missing required argument. Check ?MIRAPipeline.")
    }
    
    
    bigBin = lapply(X=BSDTList, FUN=aggregateMethyl, GRList=regionGRList, 
                    binNum=binNum, minBaseCovPerBin=minBaseCovPerBin)
    bigBinDT1 = rbindNamedList(bigBin)
    setkey(annoDT, sampleName)
    setkey(bigBinDT1, sampleName)
    bigBinDT1 = merge(bigBinDT1, annoDT, all.x=TRUE)
    profilePlot = plotMIRAProfiles(bigBinDT1) # unnormalized
    # save(bigBinDT1, file="TFbins_healthy_enhancers.RData")
    
    bigBinDT1_raw = copy(bigBinDT1)
    bigBinDT1[, methylProp := methylProp - min(methylProp) + methylNormalizationVal, 
              by=.(featureID, sampleName)]
    profilePlot2 = plotMIRAProfiles(bigBinDT1) # normalized version
    sampleScores = bigBinDT1[, .(score = calcMIRAScore(methylProp, 
                                                       shoulderShift="auto")), 
                             by=.(featureID, sampleName)]
    setkey(annoDT, sampleName)
    setkey(sampleScores, sampleName)
    sampleScores = merge(sampleScores, annoDT, all.x=TRUE)
    scorePlot = plotMIRAScores(sampleScores)
    
    resultsList = list(bigBinDT1_raw, 
                       sampleScores, 
                       plotsList=list(profilePlot, profilePlot2, scorePlot))
    return(resultsList)
}

# Get regions in each region set (GRList) that have any overlap with a region 
# set of interest (intGR)
#
# @param GRList a list of GRanges objects
# @param intGR A single Genomic Ranges object to be intersected with all the 
#           Genomic Ranges in GR
# @param removeOL If TRUE, remove the overlapping regions and get
#       everything else
# @return newGRList The result of the intersection of intGR with each 
# Genomic Ranges instance in GR
getOLRegions <- function(GRList, intGR, removeOL=FALSE) {
    
    # converting data.table to GRanges if appropriate
    if ("data.table" %in% class(GRList[[1]])) {
        GRList = lapply(X = GRList, MIRA:::dtToGr)
    }
    if ("data.table" %in% class(intGR)) {
        intGR = MIRA:::dtToGr(intGR)
    }
    olList = lapply(X = GRList, FUN = function(x) findOverlaps(intGR, x))
    subjectHitsList = lapply(olList, FUN = function(x) sort(unique(subjectHits(x))))
    
    # intGR
    queryHitsList = lapply(olList, FUN = function(x) queryHits(x))
    # overlap any regions in GRList
    overlapAny = sort(unique(unlist(queryHitsList)))
    
    
    if (removeOL) {
        # remove overlapping regions (get non-overlapping regions)
        newGRList = mapply(FUN = function(x, ind) x[!(1:length(x) %in% ind)], 
                           GRList, subjectHitsList)
    } else { 
        # get only overlapping regions
        newGRList = mapply(FUN = function(x, ind) x[ind], GRList, 
                           subjectHitsList)
    }
    
    outputList = list(GRList=newGRList, 
                      intGROverlapInd=overlapAny)
    
    return(outputList)
}

# do dimensionality reduction on CpGs in a given region set
# center = TRUE, scale.=FALSE
# drMethod "pca" or "tsne"
# if samplesAsDimensions=TRUE, transpose methylData before dimensional reduction
# to treat samples as the dimensions and CpGs as the observations. Still
# subset to CpGs in the region set of interest. PCA is still centered (or scaled)
# for this option.
# @param ... optional arguments to prcomp or Rtsne

dimRedOnRS = function(regionSet, 
                      methylData, 
                      mCoord, 
                      drMethod = "pca", 
                      samplesAsDimensions=FALSE, ...) {
    
    # test for appropriateness of inputs/right format
    if (is(mCoord, "GRanges")) {
        coordGR <- mCoord
    } else if (is(mCoord, "data.frame")) {
        # UPDATE: does the work on data.frames that are not data.tables?
        coordGR <- COCOA:::dtToGr(mCoord)
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    
    olList <- findOverlaps(query = regionSet, subject = coordGR)
    # regionHitInd <- sort(unique(queryHits(olList)))
    cytosineHitInd <- sort(unique(subjectHits(olList)))
    thisRSMData <- t(methylData[cytosineHitInd, ])
    
    
    if (!samplesAsDimensions) { # if samples are observations (normal PCA/tSNE) 
        if (drMethod == "pca") {
            # default is center = TRUE, scale. = FALSE
            drRes = prcomp(thisRSMData, ...)
        } else if (drMethod == "tsne") {
            if (nrow(thisRSMData) < 30) {
                # perplexity must be at least as many as the data points?
                perplexity = floor(1/3 * (nrow(thisRSMData) - 1))
                drRes = Rtsne::Rtsne(thisRSMData, perplexity = perplexity)
            } else {
                drRes = Rtsne::Rtsne(thisRSMData, ...)
            }
        } else {
            error("invalid drMethod")
        }
    } else {
        if (drMethod == "pca") {
            # default is center = TRUE, scale. = FALSE
            drRes = prcomp(t(thisRSMData), ...)
        } else if (drMethod == "tsne") {
            
            if (nrow(t(thisRSMData)) < 30) {
                # perplexity must be at least as many as the data points?
                perplexity = floor(1/3 * (nrow(t(thisRSMData)) - 1))
                drRes = Rtsne::Rtsne(t(thisRSMData), perplexity = perplexity)
            } else {
                drRes = Rtsne::Rtsne(t(thisRSMData), ...)
            }
            
            
        } else {
            error("invalid drMethod")
        }
    }
    
    return(drRes)
    
}


# adds underscores to string if it doesn't already have them
# default is before and after but you can also do only "before" or only "after"
# with the side parameter
addUnderscore <- function(String, side="both") {
    if (!is.null(String)) {
        if (side == "both") {
            side <- c("left", "right")
        }
        if ("left" %in% side) {
            if (!grepl(pattern = "^_", x = String)) {
                String <- paste0("_", String)
            }
        }
        if ("right" %in% side) {
            if (!grepl(pattern = "_$", x = String)) {
                String <- paste0(String, "_")
            }
        }
    }
    return(String)
}


##### Functions for setting/controlling environment (similar to projectInit) ##
#' return input string, appended by 
getPlotSubdir = function(filename="") {
    fullPath = paste0(Sys.getenv("PLOTS"), "/", Sys.getenv("PLOT.SUBDIR"), filename)
    return(fullPath)
}



################### Functions for working with/making metadata ###############
# includes filters based on metadata

#' store information about objects from a script, 
#' each script has a different notes file
scriptNotes = function(scriptID, objectID, value = "", description = "") {
    if (!dir.exists(ffamlProc("notes"))) {
        dir.create(ffamlProc("notes"))
    }
    outFile = ffamlProc(paste0("notes", scriptID, ".txt"))
    if (!file.exists(outFile)) {
        write(paste0("Notes on objects for the ", scriptID, " script."),
              file=outFile, sep = "\t")
        write(c("objectID", "value", "description", "date"), ncolumns = 4, append = TRUE,
              file=outFile, sep = "\t")
    }
    write(c(objectID, as.character(value), description, as.character(Sys.Date())), 
          file = outFile, 
          append = TRUE, 
          ncolumns = 4,
          sep = "\t")
}

# helper function to make log files
updateLog = function(object, logFile=getOption("LOGFILE")) {
    # initialize log file if it doesn't exist
    if (!file.exists(logFile)) {
        write.table(data.table::data.table("object", "nrow", "ncol", "class"), 
                    file=logFile, row.names = FALSE, 
                    col.names = FALSE, quote = FALSE, sep = "\t")
    }
    write.table(data.frame(deparse(substitute(object)), 
                           as.character(nrow(object)), 
                           as.character(ncol(object)), 
                           as.character(class(object))),
                file=logFile, row.names = FALSE, 
                col.names = FALSE, quote = FALSE, sep = "\t", append=TRUE)
    
}



# 
# metadata should have a description column describing the source of the data
filterFetal = function(metadata) {
    fInd = grep(pattern="fetal", x = metadata$description, ignore.case=TRUE)
    fInd2 = grep(pattern="foetal", x = metadata$description, ignore.case=TRUE)
    fInd = sort(unique(c(fInd, fInd2)), decreasing = FALSE)
    return(fInd)
}


################### Functions for working with/making region sets ##############

#' @param regionSets a list of data.table regions. Each data.table should have only
#' three columns: chr, start, end in that order although those names are not required
#' returns a named list of resized Genomic Ranges objects 
prepMIRARegions = function(regionSets, newSize = 4000, rsNames) {
    lapply(regionSets, function(x) setnames(x, c("chr", "start", "end")))
    regionSets = lapply(regionSets, MIRA:::dtToGR)
    regionSets = lapply(FUN = function(x) GenomicRanges::resize(x, width = newSize, fix="center"), X = regionSets)
    names(regionSets) = rsNames
    return(regionSets)
}

# takes cpgIsland region GRL
# NOTE: shores for that same cpg island are not together in sequence of the regions
# but this should not matter if using for aggregation
#' @param regionType "shore" or "shelf"
getShoreShelf = function(regionsGR, regionType = "shore") {
    if (regionType == "shore") {
        before = IRanges(start(regionsGR) - 2001, start(regionsGR) - 1)
        beforeGR = GRanges(seqnames = seqnames(regionsGR), ranges = before)
        
        after = IRanges(end(regionsGR) + 1, end(regionsGR) + 2001)
        afterGR = GRanges(seqnames=seqnames(regionsGR), ranges = after)
        newGR = c(beforeGR, afterGR)
    } else {
        # cpg shelves
        before = IRanges(start(regionsGR) - 4002, start(regionsGR) - 2002)
        beforeGR = GRanges(seqnames = seqnames(regionsGR), ranges = before)
        
        after = IRanges(end(regionsGR) + 2002, end(regionsGR) + 4002)
        afterGR = GRanges(seqnames=seqnames(regionsGR), ranges = after)
        newGR = c(beforeGR, afterGR)        
        
    }
}

############################## Visualization functions #########################

##################A function to color cluster plot by metadata variables
# see whether samples are clustering together based on a certain variable
# compare whether two variables cooccur by comparing differently colored plots

#' @param clusteredDF A dataframe with where samples are rows, 
#' features/attributes are columns.
#' This object should have the values for samples from clustering algorithms
#' (eg PCA components) and have been merged by sample ID with the metadata table 
#' for those samples.
#' First two columns are the default columns to plot
#' on the x and y axis. Other columns should be metadata from the samples
#' which can be used to assign colors.
#' @param plotCols The names of the columns to plot on the x and y axis.
#' By default, plot first two columns.
#' A character vector.
#' @param colorByCols The names of the columns to color the plot by. 
#' There will be a separate plot for each coloring so total number 
#' of plots = length(colorByCols). A character vector.
#' @param chosenPalette only can use palettes from ColorBrewer but the code
#' could be changed to allow more, default "Set1"
colorClusterPlots = function(clusteredDF, plotCols=colnames(clusteredDF)[1:2], colorByCols, alphaVal="auto") {
    
    if (is(clusteredDF, "matrix")) {
        clusteredDF = as.data.frame(clusteredDF)
    }
    
    # making sure the chosen columns are valid
    absentXYCols = !(plotCols %in% colnames(clusteredDF))
    if (any(absentXYCols)) {
        warning("One or both of given column names for plotCols were not present. Using first two columns for plotCols instead.")
        plotCols=colnames(clusteredDF)[1:2]
    }
    
    # setting alpha if there are many points
    if (alphaVal == "auto") {
        sampleNum = nrow(clusteredDF)
        if (sampleNum < 30) {
            alphaVal = 1
        } else {
            # the higher the numerator, the slower the decrease with sampleNum
            # 70 is there to shift numbers so that for 30 samples alphaVal = 1
            alphaVal = round(200 / (sampleNum + 170), 2)
        }
        
    }
    
    
    # columns to color by
    absentCols = !(colorByCols %in% colnames(clusteredDF))
    # if any are absent
    if (any(absentCols)) {
        warning("Some of the given column names were not present in colorByCols. These will be skipped.")
        # getting rid of absent columns
        colorByCols = colorByCols[!absentCols]
    }
    
    # setting up the plot
    basePlot = ggplot(data = clusteredDF, mapping = aes_string(x=plotCols[1], y=plotCols[2])) # + 
    # scale_color_brewer(palette = chosenPalette)
    
    plotList = list()
    for (i in seq_along(colorByCols)) {
        plotList[[i]] = basePlot + geom_point(aes_string(color=colorByCols[i]),alpha=alphaVal) +
            theme(aspect.ratio = 1)
    }
    
    
    # two columns, number of rows depends on how many variables to color by
    multiColPlots = marrangeGrob(plotList, ncol = 2, nrow = 2)
    return(multiColPlots)
}


# pairwise PC plots colored by colorByCols, all in one pdf if only one colorByCols
# otherwise, separate pdf for each pc combination
plotPairwiseColPCs <- function(pcaWithAnno, pcsToPlot=1:6, colorByCols, plotDir=getwd(), nameString=NULL,
                               alphaVal="auto") {
    
    if(!dir.exists(plotDir)) {
        dir.create(plotDir, recursive=TRUE)
    }
    
    
    if (!is.null(nameString)) {
        nameString <- addUnderscore(nameString)
    }
    pcComb = combn(pcsToPlot, 2)
    
    # if multiple colorByCols, plot each PC combination in a separate file
    # if only one colorbyCols, plot all in one file
    if (length(colorByCols) == 1) {
        plotCols = paste0("PC", pcsToPlot)
        # making sure the chosen columns are valid
        absentXYCols = !(plotCols %in% colnames(pcaWithAnno))
        if (any(absentXYCols)) {
            warning("One or both of given column names for plotCols were not present. Using first two columns for plotCols instead.")
            plotCols=colnames(pcaWithAnno)[1:2]
        }
        
        # setting alpha if there are many points
        if (alphaVal == "auto") {
            sampleNum = nrow(pcaWithAnno)
            if (sampleNum < 30) {
                alphaVal = 1
            } else {
                # the higher the numerator, the slower the decrease with sampleNum
                # 70 is there to shift numbers so that for 30 samples alphaVal = 1
                alphaVal = round(200 / (sampleNum + 170), 2)
            }
            
        }
        
        
        # columns to color by
        absentCols = !(colorByCols %in% colnames(pcaWithAnno))
        # if any are absent
        if (any(absentCols)) {
            warning("Some of the given column names were not present in colorByCols. These will be skipped.")
            # getting rid of absent columns
            colorByCols = colorByCols[!absentCols]
        }
        
        # setting up the plot
        basePlot =  # + 
        # scale_color_brewer(palette = chosenPalette)
        
        plotList = list()
        for (i in 1:ncol(pcComb)) {
            
            plotpc1 = pcComb[ 1, i]
            plotpc2 = pcComb[ 2, i]
            plotList[[i]] = ggplot(data = pcaWithAnno, mapping = aes_string(x=paste0("PC", plotpc1), y=paste0("PC", plotpc2))) + geom_point(aes_string(color=colorByCols[1]),alpha=alphaVal) +
                theme(aspect.ratio = 1)
        }
        multiColPlots = marrangeGrob(plotList, ncol = 2, nrow = 2)
        ggplot2::ggsave(filename=paste0(plotDir, "/multiColorPCAPlots_allPCs",
                                        nameString, 
                                        ".pdf"), plot = multiColPlots, device = "pdf",
                        limitsize=FALSE)
        
        
    } else if (length(colorByCols) > 1) {

        # do by PC combination instead of by colorByCols
        # # once for each colorByCols
        # for (j in seq_along(colorByCols)) {
        # 
        
            # once for each pcCombination
            for (i in 1:ncol(pcComb)) {
                
                plotpc1 = pcComb[ 1, i]
                plotpc2 = pcComb[ 2, i]
                multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                                       plotCols = c(paste0("PC", plotpc1), paste0("PC", plotpc2)), 
                                                       colorByCols=colorByCols)
                ggplot2::ggsave(filename=paste0(plotDir, "/multiColorPCAPlots_",
                                                plotpc1, "x", plotpc2, nameString,
                                                ".pdf"), plot = multiColorPCAPlots, device = "pdf",
                                limitsize=FALSE)
            }
        #}
    } else {
            warning("colorByCols had length 0.")
    }
}
    
    
    


#' pairwise PC plots with text IDs, pcaWithAnno must have subject_ID column which
#' will be used to label samples on plot
plotPairwiseTxtPCs <- function(pcaWithAnno, pcsToPlot=1:6, plotDir=getwd(), nameString=NULL) {
    
    if(!dir.exists(plotDir)) {
        dir.create(plotDir, recursive = TRUE)
    }
    
    
    if (!is.null(nameString)) {
        nameString <- addUnderscore(nameString)
    }
    pcComb = combn(pcsToPlot, 2)
    
    # once for each pcCombination
    for (i in 1:ncol(pcComb)) {
        
        plotpc1 = pcComb[ 1, i]
        plotpc1 = paste0("PC", plotpc1)
        plotpc2 = pcComb[ 2, i]
        plotpc2 = paste0("PC", plotpc2)
        
        
        ggPCATextPlot = ggplot(data = pcaWithAnno, mapping = aes_string(x=plotpc1, y=plotpc2)) + 
            geom_text(label=pcaWithAnno$subject_ID) + theme(aspect.ratio=1)
        ggplot2::ggsave(filename=paste0(plotDir, "/pcaTextPlot_",
                                        plotpc1, "x", plotpc2, nameString, 
                                        ".pdf"), plot = ggPCATextPlot, device = "pdf",
                        limitsize=FALSE)
    }
    
}

percentOL = function(query, subject) {
    OL = findOverlaps(query=query, subject=subject)
    return(length(unique(subjectHits(OL))) / length(subject))
}

################################################################################
# add BSAggregate from RGenomeUtils so I won't have to use that package

#' BSaggregate -- Aggregate a BSDT across regions or region groups,
#' for multiple samples at a time.
#' This function is as BScombineByRegion, but can handle not only multiple
#' samples in BSDT, but also simultaneously multiple region sets by passing
#' a regionsGRL (GRangesList object).
#' you can use jExpr to do other functions.

#' Given a bisulfite data table as input, with an identifier column for
#' different samples; plus a GRanges objects with regions to aggregate.
#'
#' @param BSDT The bisulfite data.table (output from one of the parsing
#' functions for methylation calls) that you wish to aggregate. It can
#' be a combined table, with individual samples identified by column passed
#' to splitFactor.
#' @param regionsGRL Regions across which you want to aggregate.
#' @param excludeGR A GenomicRanges object with regions you want to 
#' exclude from the aggregation function. These regions will be eliminated
#' from the input table and not counted.
#' @param jExpr You can pass a custom command in the j slot to data.table
#' specifying which columns to aggregate, and which functions to use. You
#' can use buildJ() to build a jExpr argument easily.
#' @param byRegionGroup You can aggregate by regionID or by regionGroupID; 
#' this reflects the regionsGRL that you pass; by default, BSAggregate will
#' aggregate each region individually -- scores will then be contiguous, and
#' the output is 1 row per region.
#' Turn on this flag to aggregate across all region groups, making the result
#' uncontiguous, and resulting in 1 row per *region group*.
#'
#' @export
BSAggregate = function(BSDT, regionsGRL, excludeGR=NULL, regionsGRL.length = NULL, splitFactor=NULL, keepCols=NULL, 
                       sumCols=NULL, jExpr=NULL, byRegionGroup=FALSE, keep.na=FALSE) {
    
    # Assert that regionsGRL is a GRL.
    # If regionsGRL is given as a GRanges, we convert to GRL
    if( "GRanges" %in% class(regionsGRL)) {
        regionsGRL = GRangesList(regionsGRL)
    } else if (! "GRangesList" %in% class(regionsGRL)) {
        stop("regionsGRL is not a GRanges or GRangesList object")
    }
    
    if(! is.null(excludeGR)) {
        BSDT = BSFilter(BSDT, minReads=0, excludeGR)
    }
    
    bsgr = BSdtToGRanges(list(BSDT))
    
    additionalColNames = setdiff(colnames(BSDT), c("chr","start", "end","hitCount","readCount", splitFactor))
    
    colModes = sapply(BSDT,mode)
    if (is.null(sumCols)) {
        sumCols = setdiff(colnames(BSDT),c("chr", "start", "end", "strand", splitFactor, keepCols))
        # Restrict to numeric columns.		
        sumCols = intersect(sumCols, names(colModes[which(colModes == "numeric")]))
        
    }
    # It's required to do a findoverlaps on each region individually,
    # Not on a GRL, because of the way overlaps with GRLs work. So,
    # we must convert the GRL to a GR, but we must keep track of which
    # regions came from which group.
    regionsGR = unlist(regionsGRL)
    
    if(is.null(regionsGRL.length)) {
        if (length(regionsGRL) > 100) {
            message("BSAggregate: Calculating sizes. You can speed this up by supplying a regionsGRL.length vector...", appendLF=FALSE)
        }
        regionsGRL.length = sapply(regionsGRL, length)
        message("Done counting regionsGRL lengths.")
    }
    
    # Build a table to keep track of which regions belong to which group
    region2group = data.table(
        regionID=1:length(regionsGR), 
        chr=as.vector(seqnames(regionsGR)), 
        start=as.vector(start(regionsGR)), 
        end=as.vector(end(regionsGR)),
        withinGroupID= as.vector(unlist(sapply(regionsGRL.length, seq))),
        regionGroupID=rep(1:length(regionsGRL), regionsGRL.length))
    setkey(region2group, regionID)
    
    
    message("Finding overlaps...")
    fo = findOverlaps(bsgr[[1]], regionsGR)
    
    setkey(BSDT, chr, start)
    # Gut check:
    # stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))
    
    message("Setting regionIDs...")
    BSDT = BSDT[queryHits(fo),] #restrict the table to CpGs in any region.
    
    if (NROW(BSDT) < 1) {
        warning("No BSDT sites in the given region list; please expand your regionsGRL")
        return(NULL)
    }
    
    BSDT[,regionID:=subjectHits(fo)] #record which region they overlapped.
    #BSDT[queryHits(fo),regionID:=subjectHits(fo)]
    #if (!keep.na) {
    #	BSDT = BSDT[queryHits(fo),]
    #}
    
    if (is.null(jExpr)) {
        cols=c(sumCols, keepCols)
        funcs = c(rep("sum", length(sumCols)), rep("unique", length(keepCols)))
        jExpr = buildJ(cols, funcs)
    }
    message("jExpr: ", jExpr)
    
    # Define aggregation column. aggregate by region or by region group?
    if (byRegionGroup) {
        agCol = "regionGroupID"
    } else {
        agCol = "regionID" # Default
    }
    
    # Build the by string
    if (is.null(splitFactor)) {
        byString = paste0("list(regionID)")
    } else {
        byString = paste0("list(", paste("regionID", paste0(splitFactor, ""), collapse=", ", sep=", "), ")")
    }
    
    # Now actually do the aggregate:
    message("Combining...")
    bsCombined = BSDT[,eval(parse(text=jExpr)), by=eval(parse(text=byString))]
    setkey(bsCombined, regionID)
    # Now aggregate across groups.
    # I do this in 2 steps to avoid assigning regions to groups,
    # which takes awhile. I think this preserve memory and is faster.
    
    # Define aggregation column. aggregate by region or by region group?
    if (byRegionGroup) {
        # must set allow=TRUE here in case there are multiple IDs (splitCol)
        bsCombined[region2group, regionGroupID:=regionGroupID, allow=TRUE]
        if (! is.null(splitFactor) ) { 
            byStringGroup = paste0("list(", paste("regionGroupID", paste0(splitFactor, collapse=", "), sep=", "), ")")
        } else {
            byStringGroup = "list(regionGroupID)"
        }
        bsCombined=bsCombined[,eval(parse(text=jExpr)), by=eval(parse(text=byStringGroup))]
        return(bsCombined)
    } else {
        e = region2group[bsCombined,]
        setkey(e, regionID)
        return(e)
    }
    # WARNING: There are now 2^2 ways to aggregate, sum vs mean
    # at each level: across regions, then across region sets. THis
    # doesn't give you a choice at this point. 
}
    
################################################################################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(plots, cols=1, layout=NULL) {
    library(grid)
    
    
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

