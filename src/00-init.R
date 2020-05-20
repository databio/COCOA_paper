# load libraries and prepare environment for other scripts

# for rivanna only
.libPaths(c("/home/jtl2hk/containerFiles/3.6/", .libPaths()))
.libPaths(c(.libPaths(), "/home/jtl2hk/R/x86_64-pc-linux-gnu-library/3.6"))

require(dplyr)
require(tidyr)
require(LOLA)
require(simpleCache)
require(data.table)
require(ggplot2)
require(GenomicRanges) # GRangesList, resize
# require(caret)
# require(RGenomeUtils)
require(MIRA)
require(ComplexHeatmap)
require(gridExtra) #marrangeGrob for colorClusterPlots()
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-genericFunctions.R" )) 
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-dataProcessingFunctions.R" ))
require(MultiAssayExperiment)
require(folderfun)
require(COCOA)



# options(stringsAsFactors = FALSE)

# source(paste0(Sys.getenv("CODE"), "COCOA/R/COCOA.R"))
# source(paste0(Sys.getenv("CODE"), "COCOA/R/visualization.R"))
# devtools::load_all(ffCode("COCOA/"))


brcaMetadata = fread(paste0(Sys.getenv("CODE"), 
                         "COCOA_paper/metadata/brca_metadata.csv"))
# only keep patients who have definitive status for ER and PGR
brcaMetadata = brcaMetadata[brcaMetadata$ER_status %in% 
                                                     c("Positive", "Negative"), ]
brcaMetadata = brcaMetadata[brcaMetadata$PGR_status %in% 
                                c("Positive", "Negative"), ]


# indexed clinical data (brca_clinical_metadata.tsv) has more up to date follow up info
brcaMetadata2 = read.table(file = paste0(Sys.getenv("CODE"), 
                                         "COCOA_paper/metadata/brca_clinical_metadata.tsv"), 
                           sep = "\t", header = TRUE)

patientMetadata = merge(brcaMetadata, brcaMetadata2[, c("bcr_patient_barcode", 
                                                       "vital_status", 
                                                       "days_to_death", 
                                                       "days_to_last_follow_up")],
                        by.x="subject_ID", 
                        by.y="bcr_patient_barcode", all.x=TRUE)

row.names(patientMetadata) <- patientMetadata$subject_ID




setSubPlotDir = function(dirName) {
    Sys.setenv("SUBPLOTDIR"=paste0(dirName, "/"))
}

dirPlot = function(plotFile) {
    if (!is.na(Sys.getenv("SUBPLOTDIR", unset = NA))) {
        return(paste0(Sys.getenv("PLOTS"), Sys.getenv("SUBPLOTDIR"), plotFile))    
    } else(
        return(paste0(Sys.getenv("PLOTS"), plotFile))
    )
    
}

# the AML init script will set a different plots directory
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
setff("Plot", Sys.getenv("PLOTS"))
setff("Proc", Sys.getenv("PROCESSED"))
setff("Code", paste0(Sys.getenv("CODE")))
setff("ProjCode", paste0(Sys.getenv("CODE"), "COCOA_paper/src/"))
setff("Data", Sys.getenv("DATA"))
setff("Sheets", paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/"))

dirCode = function(.file="") {
    return(paste0(Sys.getenv("CODE"), "COCOA_paper/", .file))
}


createPlotSubdir <- function(plotSubdir) {
    if (!dir.exists(ffPlot(plotSubdir))) {
        dir.create(ffPlot(plotSubdir))
    }
}

# for ggplot2
theme_set(theme_classic() + 
              theme(text = element_text(colour = "black", size = 12.5),
                    axis.text=element_text(colour = "black", size = 12.5),
                    axis.ticks = element_line(colour = "black")))
#theme(axis.text = element_text(colour = "black", size = 15), axis.ticks = element_line(colour = "black")))


# set environment
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))

#######################################################################
# function to add rank col

addRankCol = function(dataDF, colToRank, newColName="rank", decreasing=FALSE) {
    
    origColNames = colnames(dataDF)
    origInd = 1:nrow(dataDF)
    indDF = data.frame(origInd = origInd[order(dataDF[, colToRank], 
                                               decreasing = decreasing)])
    indDF$rank = 1:nrow(indDF)
    indDF = indDF[order(indDF$origInd, decreasing=FALSE), ]
    dataDF = cbind(dataDF, indDF$rank)
    colnames(dataDF) = c(origColNames, newColName)
    # updated DF
    return(dataDF)
}

# rows in rsScores should correspond to rsCollection
screenOutRoadmap = function(rsScores, rsCollection, patternSearchCol="rsName", keepPattern=NULL) {
    if (!is.null(keepPattern)) {
        keepInd =  (rsCollection != "roadmap_epigenomics") | (grepl(pattern = keepPattern, x = rsScores[, patternSearchCol], ignore.case = TRUE))
    } else {
        keepInd = (rsCollection != "roadmap_epigenomics")
    }
    
    return(keepInd)
}

conAssign <- function(vName, value, .env=currentEnv) {
    currentEnv = parent.frame(n=1) 
    
    # pull in outside env to check for existence
    if (!exists(vName, envir = .env)) {
        # assign to parent env
        assign(x = vName, value = value, envir = .env)
    }
}

############ functions to add to COCOA ##############################
# ggplot version of rs concentration
# 1 row per region set, column for rank in a given PC, 0/1 column for ER or not

# @param useGlobalTheme logical. If global ggplot theme has been set, it will
# be used if useGlobalTheme=TRUE. plot + theme_get() 
# TODO: add check that scoreColNames are present as columns of rsScores
plotRSConcentration <- function(rsScores, scoreColName="PC1", 
                                colsToSearch = c("rsName", "rsDescription"), 
                                pattern, percent = FALSE, 
                                binwidth=50, useGlobalTheme=TRUE) {
    # breaks
    rsScores = as.data.frame(rsScores)
    rsRankInd = rsRankingIndex(rsScores=rsScores, signalCol=scoreColName)
    
    
    rsInd = rep(FALSE, nrow(rsScores))
    for (i in seq_along(colsToSearch)) {
        rsInd = rsInd | grepl(pattern = pattern, x = rsScores[, colsToSearch[i]], ignore.case = TRUE)
    }
    
    rsScores$ofInterest = rsInd
    # reorder logical vector for each PC
    ofInterestDF = as.data.frame(mapply(FUN = function(ind) rsInd[as.matrix(rsRankInd)[, ind]], 
                          ind=seq_along(scoreColName))) 
    colnames(ofInterestDF) <- colnames(rsRankInd)
    ofInterestDF$rsRank = 1:nrow(ofInterestDF)
    # reshaping to long format
    ofInterestDF = tidyr::gather(ofInterestDF, key="PC", value="of_interest", scoreColName)
    if (sum(scoreColName %in% colnames(rsScores)) == 1) {
        categoryDistPlot = ggplot(ofInterestDF, aes(x=rsRank, weight=of_interest)) + 
            geom_histogram(binwidth=binwidth, boundary=0, closed="right") + theme_classic() + xlab("Region set rank") +
            ylab(paste0("Number of region sets (binwidth=", binwidth, ")")) + ggtitle(scoreColName[scoreColName %in% colnames(rsScores)])
    } else {
        categoryDistPlot = ggplot(ofInterestDF, aes(x=rsRank, weight=of_interest)) + 
            geom_histogram(binwidth=binwidth, boundary=0, closed="right") + theme_classic() + xlab("Region set rank") +
            ylab(paste0("Number of region sets (binwidth=", binwidth, ")")) + facet_wrap(~PC)
    }
    
    if (useGlobalTheme) {
        categoryDistPlot = categoryDistPlot + theme_get()    
    }
    
    
    return(categoryDistPlot)
}

#################################################################################

# can order region sets by raw score or p value
# rawScores and pVals should have same column names: colsToAnnotate
# @param pVals data.frame
# @param rankBy character. "pVals" ranks by pVal and then rawScore

#' @return 
formattedCOCOAScores  <- function(rawScores,
                                  colsToAnnotate=paste0("PC", 1:10), 
                                  numTopRS=50, pVals=NULL, rankBy=c("rawScores", "pVals")) {
    
    topRSN = numTopRS # this many top RS for each colsToAnnotate
    pRankedScores = rawScores
    
    # combine
    if (!is.null(pVals)) {
        if (pVals$rsName != rawScores$rsName) {
            stop("rawScores and pVals must be in the same order and same length.")
        }
        if (is.data.table(pVals)) {
            pVals = as.data.frame(pVals)
        }
        
        pVals = pVals[, colsToAnnotate]
        colnames(pVals) = paste0(colnames(pVals), "_pValue")
        
        pRankedScores = cbind(pRankedScores, pVals)
    }
    
    pRankedScores$index = 1:nrow(pRankedScores)
    
    # get top region sets for each colsToAnnotate based on "rankBy"
    topRSZAnnoList = list()
    
    for (i in seq_along(colsToAnnotate)) {
        
        if (is.null(pVals)) {
            
            theseTopInd = dplyr::arrange(pRankedScores,
                                         desc(get(colsToAnnotate[i])))$index[1:topRSN]
            # thesePValRanks = order(pRankedScores[, paste0(colsToAnnotate[i])], decreasing = FALSE)
            # pRankedScores$index[thesePValRanks]
            
            topRSZAnnoList[[i]] = data.frame(pRankedScores[theseTopInd, c("rsName", "rsDescription", colsToAnnotate[i],
                                                                          "signalCoverage", "regionSetCoverage",
                                                                          "totalRegionNumber", "meanRegionSize")])
            
            names(topRSZAnnoList[[i]]) <- paste0(colsToAnnotate[i], "_", c("rsName", "rsDescription", "rsScore",
                                                                           "signalCoverage", "regionSetCoverage",
                                                                           "totalRegionNumber", "meanRegionSize"))
        } else {
            if (rankBy == "pVals") {
                theseTopInd = dplyr::arrange(pRankedScores, 
                                             desc(get(paste0(colsToAnnotate[i], "_PValue"))), 
                                             desc(get(colsToAnnotate[i])))$index[1:topRSN]
            } else {
                theseTopInd = dplyr::arrange(pRankedScores, 
                                             desc(get(colsToAnnotate[i])))$index[1:topRSN]
            }

            topRSZAnnoList[[i]] = data.frame(pRankedScores[theseTopInd, c("rsName", "rsDescription", colsToAnnotate[i], 
                                                                          paste0(colsToAnnotate[i], "_PValue"),
                                                                          "signalCoverage", "regionSetCoverage",
                                                                          "totalRegionNumber", "meanRegionSize")])
            
            names(topRSZAnnoList[[i]]) <- paste0(colsToAnnotate[i], "_", c("rsName", "rsDescription", "rsScore",
                                                                           "PValue",
                                                                           "signalCoverage", "regionSetCoverage",
                                                                           "totalRegionNumber", "meanRegionSize"))
        }
    }
    
    
    return(topRSZAnnoList)
    
    # # p val cutoffs
    # sigCutoff = 0.05
    # trendCutoff = 0.1
    # 
    # # sort region sets according to p-value/z-score groups (significant, trending 
    # # toward significant, nonsignificant), then by average correlation score
    # pValGroups = as.data.frame(gPValDF)[, colsToAnnotate, drop=FALSE]
    # 
    # sigInd = pValGroups <= sigCutoff
    # notSigInd = pValGroups > trendCutoff
    # trendInd = (pValGroups <= trendCutoff) & (pValGroups > sigCutoff)
    # 
    # pValGroups[sigInd] = 1
    # pValGroups[notSigInd] = -1
    # pValGroups[trendInd] = 0
    # 
    # 
    # pValGroups = cbind(pValGroups, gPValDF[, c("rsName", "rsDescription")])
    # colnames(pValGroups) = paste0(colnames(pValGroups), "_PValGroup")
    # 
    # pRankedScores = cbind(realRSScores, pValGroups)
    # # View(dplyr::arrange(pRankedScores, desc(LF4_Z), desc(LF4)))
    # pRankedScores$index = 1:nrow(pRankedScores)
    # gPValDF2 = as.data.frame(gPValDF)[, colsToAnnotate, drop=FALSE]
    # colnames(gPValDF2) <- paste0(colnames(gPValDF2), "_PVal")
    # pRankedScores = cbind(pRankedScores, gPValDF2)
    # realRSScores = rsScores
}


# @param rsScore numeric A single number that is the score of a given region set
# for this PC ("pc")
# @param list Each item in this list is a matrix that has 
# @param numeric The number of regions from the region set that overlapped 
# with any of the data ("genomicSignal")
# returns a list of two matrices, one matrix for upper bounds 
# and one matrix for lower bounds
getPValRange <- function(rsScore, nullDistList, pc, regionCoverage) {
    sapply(X = nullDistList, FUN = function(x) sum(x[ , pc] > rsScore) / nrow(x))
}
getPValRangeVec <- function(rsScore, nullDistList, pc, regionCoverage) {
    getPValRange(rsScore, nullDistList, pc, regionCoverage)
}



getUpperBound <- function(rsScore, nullDistList, sampleSize, pc, regionCoverage) {
    
    getUpperBoundSingle <- function(rsScore, nullDistList, 
                                    sampleSize, pc, regionCoverage) {
        belowInd = which(sampleSize < regionCoverage)
        
        # check for empty vector (none were below)
        if (length(belowInd) == 0) {
            # worst possible p value
            bound = 1
            return(bound)
        }
        
        lastBelow = belowInd[length(belowInd)]
        
        bound = sum(nullDistList[[lastBelow]][ , pc] > rsScore) / nrow(nullDistList[[lastBelow]])
        
        # if zero, upper bound is lowest possible known p val
        if (bound == 0) {
            bound = 1 / nrow(nullDistList[[lastBelow]])
        }
        
        return(bound)
    }
    
    # get upper bound for all scores in vector (rsScore and regionCoverage are vectors
    # of the same length)
    bound = mapply(FUN = function(x, y) getUpperBoundSingle(rsScore = x, 
                                                            nullDistList = nullDistList, sampleSize=sampleSize, 
                                                            pc=pc, regionCoverage = y), 
                   x = rsScore, y = regionCoverage)
    return(bound)
    
} 
getLowerBound <- function(rsScore, nullDistList, sampleSize, pc, regionCoverage) {
    
    getLowerBoundSingle <- function(rsScore, nullDistList, sampleSize, pc, regionCoverage) {
        aboveInd = which(sampleSize > regionCoverage)
        
        # check for empty vector (none were below)
        if (length(aboveInd) == 0) {
            # best possible p value
            bound = 0
            return(bound)
        }
        
        firstAbove = aboveInd[1]
        bound = sum(nullDistList[[firstAbove]][ , pc] > rsScore) / nrow(nullDistList[[firstAbove]])
        return(bound)
    }
    
    # get lower bound for all scores in vector (rsScore and regionCoverage are vectors
    # of the same length)
    bound = mapply(FUN = function(x, y) getLowerBoundSingle(rsScore = x, 
                                                            nullDistList = nullDistList, sampleSize=sampleSize, 
                                                            pc=pc, regionCoverage = y), 
                   x = rsScore, y = regionCoverage)
    return(bound)
    
} 



#####################################################################
# functions for permutations (moved to COCOA package)


################## not sure whether to add this one to COCOA
# convenience function to quickly make meta-region loading profiles
# first calculates profiles, then normalizes and plots
# output is list with 3 items:  c("grob", "metaRegionData", "plotList")
# grob: an marrangeGrob of all plots and should be saved with ggsave
# metaRegionData: unnormalized output of getMetaRegionProfile
# plotList: the plot for each region set individually (all signalCol in one plot)
#' @param aggrMethod character. A character object with the aggregation method.
#' Similar to `scoringMetric`.
#' There are different aggregation methods available for 
#' signalCoordType="singleBase" vs  signalCoordType="multiBase".
#' For "singleBase", the available scoring methods are "regionMean" and
#' "simpleMean". The default method is "regionMean".
#' For "multiBase", the scoring methods are "proportionWeightedMean" and 
#' "simpleMean". The default is "proportionWeightedMean".
#' "regionMean" is a weighted
#' average of the signal, weighted by region (absolute value of signal 
#' if absVal=TRUE). First the signal is
#' averaged within each regionSet region, 
#' then all the regions are averaged. With
#' "regionMean" score, be cautious in interpretation for
#' region sets with low number of regions that overlap signalCoord. 
#' The "simpleMean"
#' method is just the unweighted average of all (absolute) signal values that
#' overlap the given region set. For multiBase data, this includes
#' signal regions that overlap a regionSet region at all (1 base
#' overlap or more) and the signal for each overlapping region is
#' given the same weight for the average regardless of how much it overlaps. 
#' "proportionWeightedMean" is a weighted average of all signalCoord 
#' regions that overlap with regionSet regions. For each signalCoord region
#' that overlaps with a regionSet region, we calculate what proportion
#' of the regionSet region is covered. Then this proportion is used to
#' weight the signal value when calculating the mean. 
#' The denominator of the mean
#' is the sum of all the proportion overlaps. 
#' @return returns a named list of three elements: "grob", "metaRegionData", 
#' and "plotList". 1st element is the grob of all plots, which can be
#' plotted with ggplot. 2nd element is normalized output of getMetaRegionProfile.
#' 3rd element is a list where each item is a single plot.

makeMetaRegionPlots <- function(signal, signalCoord, GRList, rsNames, 
                                signalCol, binNum, 
                                # returnNormalizedVals=TRUE, 
                                aggrMethod="default", absVal=TRUE, 
                                normMethod = c("mean", "none", "zscore")) {
    
    pcProf = lapply(X = GRList, function(x) getMetaRegionProfile(signal = signal, 
                                                              signalCoord = signalCoord, 
                                                              regionSet = x, signalCol = signalCol,
                                                              binNum = binNum, 
                                                              aggrMethod=aggrMethod,
                                                              absVal = absVal))
    
    # get rid of region sets that had insufficient overlap
    notNull = !vapply(X = pcProf, FUN = is.null, FUN.VALUE = TRUE)
    pcProf = pcProf[notNull]
    rsNames = rsNames[notNull]
    
    # emptyInd = sapply()
    # # for single columns?
    # # average loading value from each PC to normalize so PCs can be compared with each other
    # if (is.numeric(brcaLoadings[, PCsToAnnotate])) {
    #     avLoad = mean(abs(brcaLoadings[, PCsToAnnotate]))
    # } else {
    #     avLoad <- apply(X=brcaLoadings[, PCsToAnnotate], 
    #                     MARGIN=2, 
    #                     FUN=function(x) mean(abs(x)))
    # }
    
    # normalize
    # pcP = lapply(pcP, FUN = function(x) t(apply(X = x, MARGIN = 1, FUN = function(y) y - c(0, avLoad))))
    # for (i in seq_along(signalCol)) {
    #     # by reference
    #     lapply(pcP, FUN = function(x) x[, 1])
    # }
    
    pcP = normalizeMRProfile(signal=signal, signalCol=signalCol, 
                                 pcProf, rsNames = rsNames, absVal=TRUE, 
                                 normMethod = normMethod)


    # for the plot scale
    maxVal = max(sapply(pcP, FUN = function(x) max(x[, loading_value])))
    minVal = min(sapply(pcP, FUN = function(x) min(x[, loading_value])))
    
    
    # stack overflow for wrapping plot title
    wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n") 
    
    
    profilePList = list()
    for (i in seq_along(pcP)) {
        
        thisRS = pcP[[i]]
        
        
        profilePList[[i]] = ggplot(data = thisRS, mapping = aes(x =binID , y = loading_value)) + 
            geom_line() + ylim(c(minVal, maxVal)) + facet_wrap(facets = "PC") + 
            ggtitle(label = wrapper(rsNames[i], width=30)) + xlab("Genome around Region Set, 14 kb") + 
            ylab("Normalized Loading Value") + 
            theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
        profilePList[[i]]
        # plot(pcP[[i]]$PC1, type="l") + title(rsNames[i])
        #xLabels = xAxisForRegionPlots2()
    }
    
    multiProfileP = marrangeGrob(profilePList, ncol = 2, nrow = 2)
    names(profilePList) <- rsNames    

    metaRegionProfileInfo = list(multiProfileP, pcProf, profilePList)
    setattr(metaRegionProfileInfo, "names", c("grob", "metaRegionData", "plotList"))
    
    return(metaRegionProfileInfo)
}

# stack overflow for wrapping plot title
wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n") 


#a data.table. ouput includes columns "PC", "loading_value"
# @param negLog logical. only applies if normMethod is making a pval, -log(x,10)
# @param format character. long or wide. If "long" is not given as "format" (long is the default),
# then wide format will be given regardless of the "format" argument
normalizeMRProfile <- function(signal, signalCol, pList, rsNames, 
                               absVal=TRUE, 
                               normMethod=c("mean", "zscore", "none", "normPVal", "empPVal"),  negLog=TRUE,
                               format="long") {
    
    pcP = copy(pList)
    
    # this check must be done before converting items of list to data.table
    # otherwise NULL will be converted to a data.table
    notNull = !vapply(X = pcP, FUN = is.null, FUN.VALUE = TRUE)
    pcP = pcP[notNull]
    rsNames = rsNames[notNull]
    pcP = lapply(pcP,FUN = as.data.table)

    # average loading value from each PC to normalize so PCs can be compared with each other
    if (absVal) {
        avLoad = apply(X = signal[, signalCol], MARGIN = 2, FUN = function(x) mean(abs(x)))
        sdLoad = apply(X = signal[, signalCol], MARGIN = 2, FUN = function(x) sd(abs(x), na.rm = TRUE))
    } else {
        avLoad = apply(X = signal[, signalCol], MARGIN = 2, FUN = function(x) mean(x))
        sdLoad = apply(X = signal[, signalCol], MARGIN = 2, FUN = function(x) sd(x, na.rm = TRUE))
    }
    
    if (normMethod == "none") {
        pcP = lapply(X = pcP, function(x) x[, signalCol, with=FALSE])
    } else if (normMethod == "mean") {
        pcP = lapply(pcP, FUN = function(x) x[, mapply(FUN = function(y, z) get(y) - z, y=signalCol, z = avLoad)])
    } else if (normMethod == "zscore") {
        pcP = lapply(pcP, FUN = function(x) as.data.table(x[, mapply(FUN = function(y, z) get(y) - z, y=signalCol, z = avLoad)]))
        pcP = lapply(pcP, FUN = function(x) x[, mapply(FUN = function(y, z) get(y) / z, y=signalCol, z = sdLoad)])
    } else if (normMethod == "normPVal") {
        pcP = lapply(pcP, FUN = function(x) as.data.table(x[, mapply(FUN = function(y, z) get(y) - z, y=signalCol, z = avLoad)]))
        pcP = lapply(pcP, FUN = function(x) as.data.table(x[, mapply(FUN = function(y, z) get(y) / z, y=signalCol, z = sdLoad)]))
        # two sided p-value (for one sided, don't multiply by two)
        pcP = lapply(X = pcP, FUN = function(x) apply(X = x, FUN = function(y) (1-pnorm(abs(y))) * 2, MARGIN = 2))
        if(negLog) {
            pcP = lapply(X = pcP, FUN = function(x) apply(X = x, FUN = function(y) -log10(y), MARGIN = 2))
        }
    } 

    pcP = lapply(pcP, FUN = function(x) data.table(binID=1:nrow(x), x))
    
    if (format == "long") {
        # convert to long format for plots
        pcP = lapply(X = pcP, FUN = function(x) tidyr::gather(data = x, key = "PC", value="loading_value", signalCol))
        pcP = lapply(X = pcP, as.data.table)
        pcP = lapply(pcP, function(x) x[, PC := factor(PC, levels = signalCol)])
    } 
    normPList = pcP
    
    return(normPList)
}

###################### visualization ##########################

# @examples dataDF = data.frame(data.frame(a=c((1:4)/10, .4), b=6:10))
# niceHist(dataDF=dataDF, colToPlot="a", xLabel="myXLab", binwidth=0.01)
# make a single ggplot histogram
niceHist = function(dataDF, colToPlot, xLabel=colToPlot, binwidth=1, boundary=0, yLabel="", plotTitle = "", 
                    ggExpr=NULL) {
    p= ggplot(dataDF, aes(get(colToPlot))) + 
        geom_histogram(binwidth=binwidth, boundary=boundary) + 
        xlab(xLabel) +ylab(yLabel) + ggtitle(plotTitle)
    # geom_histogram(aes(y=..count../sum(..count..))
    # could add arbitrary expression  # e.g. plottingExpr = "theme_classic() + (etc.)"
    if (!is.null(ggExpr)) {
        eval(parse(text=paste0("p", ggExpr)))
    } else (
        p
    )
    
    # theme_classic(), ylim(0, yMax)
}

# make one or more ggplot histograms and save them to a pdf
# @param colsToPlot character. The name of the column itself. One histogram
# for each column. 
# @param xLabels character. A more complete/understandable name/title for the 
# columns which are being plotted
multiNiceHist = function(file, dataDF, colsToPlot, xLabels=colsToPlot,
                         binwidth=1, boundary=0, yLabel="", plotTitles="", ggExpr=NULL) {
    
    # make it the same length
    if (length(plotTitles) == 1) {
        plotTitles = rep(plotTitles, length(colsToPlot))
    }
    if (length(xLabels) == 1) {
        xLabels = rep(xLabels, length(colsToPlot))
    }
    
    pdf(file = file)
    
    for (i in seq_along(colsToPlot)) {
        colToPlot = colsToPlot[i]
        p = niceHist(dataDF, colToPlot=colToPlot, xLabel=xLabels[i], 
                     binwidth=binwidth, boundary=boundary, yLabel=yLabel,
                     plotTitle=plotTitles[i], ggExpr=ggExpr)  
        print(p)
    }
    
    dev.off()
    
}


# same as plotAnnoScoreDist but add an extra annotation bar on top of the plot
# this requires using cowplot to put multiple plots together
plotAnnoScoreDist2 <- function(rsScores, colsToPlot, pattern, patternName=pattern) {
    require("cowplot")
    
    rsCorP = plotAnnoScoreDist(rsScores=rsScores, colsToPlot=colsToPlot, 
                               pattern=pattern, patternName=patternName)
    ##################
    
    rsScores$rank = order(order(as.numeric(rsScores[, colsToPlot]), decreasing = TRUE))
    
    rsGroup = rep(0, nrow(rsScores))
    for (i in seq_along(pattern)) {
        thisPatternInd = grepl(pattern = pattern[i], 
                               x = rsScores$rsName, 
                               ignore.case = TRUE)
        # earlier groups will be overwritten if they are not mutually exclusive
        # (e.g. if the same region set matches two patterns, the later group
        # will be assigned)
        rsGroup[thisPatternInd] = i
    }
    rsGroup = as.factor(rsGroup)
    # necessary so region set names will be in legend
    levels(rsGroup) = c("Other", patternName)
    rsScores$Group = rsGroup 
    
    ###############
    # as.numeric(Group)-1)
    # rsCorP + coord_flip()
    polycombStatusP = ggplot(data = rsScores, mapping = aes(x=rank, y=1, col=Group)) + 
        # geom_col(aes(col=Group), alpha=0) +
        scale_color_discrete(drop = FALSE) + theme(aspect.ratio = 0.03) + scale_x_continuous(limits = c(0, nrow(rsScores)))
    
    for (i in 2:(length(pattern) + 1)) {
        polycombStatusP = polycombStatusP + geom_col(data = rsScores[as.numeric(rsScores$Group) == i, ], 
                                                     alpha=1)#, show.legend = FALSE)
    }
    
    polycombStatusP 
    
    gg_dist_g1 = polycombStatusP
    # gg_dist_g2 = rsCorP
    gg_scatter = rsCorP
    
    # gg_dist_g2 = gg_dist_g2 + coord_flip()
    
    # Remove some duplicate axes
    gg_dist_g1 = gg_dist_g1 + theme(axis.title.x=element_blank(),
                                    axis.title.y = element_blank(),
                                    axis.text=element_blank(),
                                    axis.line=element_blank(),
                                    axis.ticks=element_blank())
    #legend.position = "bottom")
    #legend.text = element_blank())
    gg_dist_g1
    
    # Modify margin c(top, right, bottom, left) to reduce the distance between plots
    #and align G1 density with the scatterplot
    gg_dist_g1 = gg_dist_g1 + theme(plot.margin = unit(c(0.5, 0, 0, 0.5), "cm"))
    gg_scatter = gg_scatter + theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"))
    #gg_dist_g2 = gg_dist_g2 + theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))
    
    # Combine all plots together and crush graph density with rel_heights
    first_col = plot_grid(gg_dist_g1, gg_scatter, ncol = 1, rel_heights = c(1, 8), align = "v", axis="lr")
    #second_col = plot_grid(NULL, gg_dist_g2, ncol = 1, rel_heights = c(1, 3))
    # perfect = plot_grid(first_col, second_col, ncol = 2, rel_widths = c(3, 1))
    p = first_col
    p
    
    return(p)
}



#' Plot region set scores and rank, annotating groups of interest
#' 
#' @param colsToPlot character. Name of the column to plot. ###update: only one can be used
#' @param pattern character. Multiple patterns can be given as character objects in a vector.
#' Regular expressions can be used (ignore.case=TRUE though)
#' @examples 
#' data(rsScores)
#' rsScores$rsName <- c("ER", "GATA3", "ER", "GATA3", "AP1") 
#' plotAnnoScoreDist(rsScores, colsToPlot="PC1", pattern="ER")
#' plotAnnoScoreDist(rsScores, colsToPlot="PC2", pattern=c("ER", "GATA3")) + geom_point(size=5)
#' @export
plotAnnoScoreDist <- function(rsScores, colsToPlot, pattern, patternName=pattern, 
                              alpha=0.5, shape=3) {
    
    rsScores$rank = order(order(as.numeric(rsScores[, colsToPlot]), decreasing = TRUE))
    
    rsGroup = rep(0, nrow(rsScores))
    for (i in seq_along(pattern)) {
        thisPatternInd = grepl(pattern = pattern[i], 
                               x = rsScores$rsName, 
                               ignore.case = TRUE)
        # earlier groups will be overwritten if they are not mutually exclusive
        # (e.g. if the same region set matches two patterns, the later group
        # will be assigned)
        rsGroup[thisPatternInd] = i
    }
    rsGroup = as.factor(rsGroup)
    # necessary so region set names will be in legend
    levels(rsGroup) = c("Other", patternName)
    rsScores$Group = rsGroup 
    
    
    rsCorP = ggplot(data = rsScores, 
                    mapping = aes(x=rank, y=get(colsToPlot), 
                                  col=Group)) + # alpha=Group
        # geom_point(alpha=0.0, shape=3) +
        ylab("Region set score") + xlab("Region set rank") +
        scale_color_discrete(drop = FALSE) + theme(aspect.ratio = 1)
    
    
    # add each group (Other and pattern) sequentially so they will be plotted on top
    # of other points ('Other' plotted first)
    
    for (i in 1:(length(pattern) + 1)) {
        rsCorP = rsCorP + geom_point(data = rsScores[as.numeric(rsScores$Group) == i, ], alpha=alpha, shape=shape) 
    }
    
    return(rsCorP)
}



###############################################################################
# @param pMeta data.frame required for making PCA plot colored by metadata
# @param colorByCols character
#### plots that will be created
## "comparePCHeatmap"
## "methylAlongPC"
## "regionQuantileByTargetVar"
## "pcFromSubset Correlation Heatmap"
## "region set Overlapping Cytosine Proportion" (rsOLCP)
## proportion of cytosines from region set that are shared with other region set
## "meta region loading profiles" (mrLP)
## annotated PC plot

# plotSubdir # directory for this analysis. A sub directory in plotSubdir will be added for the specific data input (inputID)
# inputID # since a cache is saved, this marks what the methylation data source was
# Sys.getenv("PLOTS") should be set


cocoaMultiVis <- function(sortedRSIndDF, GRList, coordinateDT, loadingMat, rsScores,
                          mPCA, PCSTOANNOTATE, methylData, plotDir, inputID, 
                          pMeta=NULL, colorByCols=NULL,
                          makeCPCH=TRUE, 
                          PCsToAnnotate_cPCH = PCSTOANNOTATE,
                          makeMAPC=TRUE, 
                          PCsToAnnotate_mAPC = PCSTOANNOTATE[1:5], topRSToPlotNum=20,
                          makeRQBPC=TRUE, 
                          topRSInd_rQBPC = unique(unlist(sortedRSIndDF[1:10, ])),
                          PCsToAnnotate_rQBPC = PCSTOANNOTATE,
                          makePCFSCH=TRUE, 
                          topRSInd_pcFSCH = unique(unlist(sortedRSIndDF[1:10, ])), 
                          PCsToAnnotate_pcFSCH=PCSTOANNOTATE,
                          makeRSOLCP=TRUE, 
                          topRSInd_rsOLCP = unique(unlist(sortedRSIndDF[1:10, ])), 
                          makeMRLP=TRUE, 
                          topRSInd_mrLP=unique(unlist(sortedRSIndDF[1:10, ])), 
                          PCsToAnnotate_mrLP = PCSTOANNOTATE,
                          makePairwisePCPlots=FALSE) {
    
    source(paste0(Sys.getenv("CODE"), "aml_e3999/src/00-genericFunctions.R"))
    source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/COCOA_extra/R/visualization.R"))
    source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/COCOA_extra/R/analysis.R"))
    require(grid)
    require(ggplot2)
    require(COCOA)
    require(ComplexHeatmap)    
    
    # place to save plots
    plotDir = paste0(plotDir, "/")
    if (!dir.exists(plotDir)) {
        dir.create(plotDir, recursive = TRUE)
    }
    
    inputID = addUnderscore(inputID, side = "left")
    
    # make new directory in specified directory since this is specific to inputID
    plotDir = paste0(plotDir, "COCOA_plots", inputID, "/")
    if (!dir.exists(plotDir)) {
        dir.create(plotDir, recursive = TRUE)
    }
    ##################################################################################
    # comparePCHeatmap
    # visualization of enrichment score results across PCs
    # see if region set high in one PC is also high in others
    
    # number of plots = length(PCsToRankBy). one plot for each
    # multiple plots in one pdf
    # TODO: filter out low coverage region sets
    # for rsScores
    if (makeCPCH) {
        tryCatch({
        # necessary
        PCsToAnnotate_cPCH
        comparePCHeatmap(rsScores=rsScores, 
                         PCsToRankBy=PCsToAnnotate_cPCH, 
                         PCsToInclude=PCsToAnnotate_cPCH,
                         fileName=paste0(plotDir, 
                                         "rsEnrichHeatmap", inputID, ".pdf"))
        },error=function(cond) {message("Error with CPCH")}, finally = {})
    }
    ##################################################################################
    # methylAlongPC
    # looking at methylation level data at individual cytosines ordered by PC 
    # only looking at regions with high average loading scores
    # still individual cytosine methylation
    
    if (makeMAPC) {
        tryCatch({
        # one pdf for each PC given.
        for (i in seq_along(PCsToAnnotate_mAPC)) {
            
            # top region sets for this PC
            rsInd = as.numeric(as.matrix(sortedRSIndDF[1:topRSToPlotNum, PCsToAnnotate_mAPC[i]])) # original index
            
            grDevices::pdf(paste0(plotDir, "regionMethylHeatmaps", PCsToAnnotate_mAPC[i], inputID, ".pdf"), width = 11, height = 8.5 * topRSToPlotNum)
            
            # heatmap
            methylAlongPC(loadingMat=loadingMat, loadingThreshold=0.95,
                          pcScores=mPCA$x,
                          signalCoord=coordinateDT,
                          genomicSignal=methylData,
                          GRList=GRList[rsInd], orderByPC=PCsToAnnotate_mAPC[i],
                          topXRegions=50,
                          cluster_columns = TRUE)
            
            # draw(Heatmap(matrix = methylData[1:1000, 1:10]))
            # plot(methylData[1:1000, 1])
            
            dev.off()
        }
        },error=function(cond) {message("Error with MAPC")}, finally = {})
    }
    ###################################################################################
    # regionQuantileByTargetVar
    # comparing loading scores/percentiles for individual regions among PCs
    # need region sets and PCA loadings
    
    if (makeRQBPC) {
        tryCatch({
        grDevices::pdf(paste0(plotDir, 
                              "regionPercentileByPC", inputID, ".pdf"), 
                       width = 11, height = 8.5 * length(topRSInd_rQBPC))
        
        ## ranking in terms of percentiles in case there were different distributions of loading scores for each PC
        # if there are too many regions, will try to cluster and cause memory error:
        # cannot allocate vector of size X Gb,
        # fix this by decreasing maxRegionsToPlot or use cluster_rows=FALSE
        multiRegionQuantileByTargetVar(signal=loadingMat, signalCoord=coordinateDT, 
                                GRList=GRList[topRSInd_rQBPC], 
                                rsNames=paste0(rsScores$rsName[topRSInd_rQBPC], " : ", rsScores$rsDescription[topRSInd_rQBPC]), 
                                signalCol=PCsToAnnotate_rQBPC, maxRegionsToPlot = 5000,
                                cluster_rows = TRUE, absVal = TRUE)
        
        dev.off()
        },error=function(cond) {message("Error with RQBPC")}, finally = {})
    }
    ##################################################################################
    # pcFromSubset Heatmap
    # seeing whether a subset of cytosines (ie a single region set) can
    # recapitulate PC score from all cytosines
    
    if (makePCFSCH) {
        tryCatch({
        .regionSetList = GRList[topRSInd_pcFSCH] 
        regionSetName = paste0(rsScores$rsName[topRSInd_pcFSCH], " : ", rsScores$rsDescription[topRSInd_pcFSCH])
        names(.regionSetList) <-  regionSetName
        subsetCorList = lapply(X = as.list(.regionSetList), FUN = function(x) pcFromSubset(regionSet = x, 
                                                                                           mPCA = mPCA, 
                                                                                           methylData = methylData, 
                                                                                           mCoord = coordinateDT, 
                                                                                           pc = PCsToAnnotate_pcFSCH,
                                                                                           returnCor = TRUE))
        # test=list()
        # for (i in seq_along(.regionSetList)) {
        #     message(i)
        #     test[[i]] = pcFromSubset(regionSet = .regionSetList[i], 
        #                              mPCA = mPCA, 
        #                              methylData = methylData, 
        #                              mCoord = coordinateDT, 
        #                              pc = PCsToAnnotate_pcFSCH,
        #                              returnCor = TRUE)
        # }
        subsetCorMat = do.call(rbind, subsetCorList)
        colnames(subsetCorMat) <- PCsToAnnotate_pcFSCH
        
        simpleCache(paste0("subsetCorMat", inputID), {subsetCorMat}, recreate=TRUE)
        
        grDevices::pdf(paste0(plotDir, "subsetCorRSbyPC", inputID, ".pdf"), width = 8.5, 11)
        
        draw(Heatmap(matrix = subsetCorMat, col = c("black", "orange"), cluster_rows = FALSE, cluster_columns = FALSE, 
                     column_title = , cell_fun = function(.j, .i, x, y, width, height, fill, mat=subsetCorMat) {
                         grid.text(sprintf("%.2f", mat[.i, .j]), x, y, gp = gpar(fontsize = 10))
                     }))
        dev.off()
        
        
        # plotting correlation between a PC and the "PC-subset" score 
        # derived from only loading values of CpGs within a certain region set 
        # do this for top few region sets (i in the loop) for each PC
        for (j in seq_along(PCsToAnnotate_pcFSCH)) {
            grDevices::pdf(paste0(plotDir, "/subsetPC_PC_Scatter_", 
                                  PCsToAnnotate_pcFSCH[j], addUnderscore(inputID, side="left"), ".pdf"))
            for (i in 1:10) {
                # returns a "PC-subset" score for each sample
                subScores = pcFromSubset(regionSet = GRList[as.numeric(as.data.frame(sortedRSIndDF)[i, PCsToAnnotate_pcFSCH[j]])], 
                                         mPCA = mPCA, 
                                         methylData = methylData, 
                                         mCoord = coordinateDT, 
                                         pc = PCsToAnnotate_pcFSCH[j],
                                         returnCor = FALSE)    
                
                plot(x = subScores, y = mPCA$x[, PCsToAnnotate_pcFSCH[j]], xlab= "Scores from only CpGs in this region set",
                     ylab = PCsToAnnotate_pcFSCH[j], main = paste0(as.data.frame(rsScores)$rsName[as.numeric(sortedRSIndDF[i, PCsToAnnotate_pcFSCH[j]])], 
                                                                   " : ", as.data.table(rsScores)$rsDescription[as.numeric(sortedRSIndDF[i, PCsToAnnotate_pcFSCH[j]])]))   
                
            }
            dev.off()
        }
        },error=function(cond) {message("Error with PCFSCH")}, finally = {})
    }
    ################################################################################
    # seeing how much overlap there is between region sets
    # based on overlap of covered cytosines, not the regions themselves

    if (makeRSOLCP) {
        tryCatch({
        # total regions in column region sets are the denominator for the proportion
        .regionSetList = GRList[topRSInd_rsOLCP] 
        pOL = percentCOverlap(mCoord = MIRA:::dtToGr(coordinateDT), 
                              GRList = .regionSetList) 
        
        grDevices::pdf(paste0(plotDir, "topRSOverlap", inputID, ".pdf"), width = 25, height = 25)
        
        draw(Heatmap(matrix = pOL[[1]], col = c("black", "yellow"), cluster_rows = FALSE, cluster_columns = FALSE, 
                     cell_fun = function(j, i, x, y, width, height, fill, mat=pOL[[1]]) {
                         grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                     }))
        
        # numbers not included on plot squares
        # Heatmap(matrix = pOL[[1]], col = c("gray14", "gold"), cluster_rows = FALSE, cluster_columns = FALSE)
        dev.off()
        },error=function(cond) {message("Error with RSOLCP")}, finally = {})
    }
    
    ####################################################################
    # "meta region loading profiles" (mrLP)
    # meta-region plots of region surrounding regions in region set
    # check whether enrichment is specific to this region set by
    # seeing if loading values have a spike in the center of these region sets
    # compared to surrounding genome 
    
    # Error in names(profilePList) <- rsNames : 
    #     'names' attribute [62] must be the same length as the vector [57]
    
    if (makeMRLP) {
        tryCatch({
        .regionSetList = GRList[topRSInd_mrLP]
        .regionSetList = lapply(.regionSetList, resize, width = 14000, fix="center")
        .rsNames = paste0(rsScores$rsName[topRSInd_mrLP], "_:_", rsScores$rsDescription[topRSInd_mrLP])
        
        # returns a list: one item is grob, one item is list of binned data tables
        mrPlotOutput = makeMetaRegionPlots(signal = loadingMat, signalCoord = coordinateDT, GRList = .regionSetList, 
                                           rsNames = .rsNames, signalCol = PCsToAnnotate_mrLP, binNum = 21, 
                                           aggrMethod = "default")
        pcProf = mrPlotOutput$metaRegionData
        
        simpleCache(paste0("pcProf14k", inputID), {
            pcProf
        }, recreate = TRUE)
        
        ggsave(filename = paste0(plotDir,
                                 "metaRegionLoadingProfiles", 
                                 inputID, ".pdf"), plot = mrPlotOutput$grob, device = "pdf", limitsize = FALSE)
        
        # check PPARG.bed, Jaspar motifs (had 18 rows instead of 21)
        },error=function(cond) {message("Error with MRLP")}, finally = {})
    }
    
    ##########################################################################
    # make PC score plot where samples are colored by metadata
    
    # TODO: add colored TSNE or UMAP
    if (!is.null(pMeta) && !is.null(colorByCols)) {
        tryCatch({
        pcScores = mPCA$x
        sharedSamples = intersect(row.names(pcScores), pMeta$sample_name)
        row.names(pMeta) = pMeta$sample_name
        colorPlot = colorClusterPlots(clusteredDF = cbind(pcScores[sharedSamples, ], pMeta[sharedSamples, ]),
                                      plotCols = PCSTOANNOTATE[1:2], colorByCols = colorByCols)
        ggplot2::ggsave(filename=paste0(plotDir, "annoPCPlot","_", PCSTOANNOTATE[1],"_", PCSTOANNOTATE[2], "_", inputID, ".pdf"),
                        plot = colorPlot, device = "pdf",
                        limitsize=FALSE)
        if (length(PCSTOANNOTATE) > 2) {
            colorPlot = colorClusterPlots(clusteredDF = cbind(pcScores[sharedSamples, ], pMeta[sharedSamples, ]),
                                          plotCols = PCSTOANNOTATE[c(1,3)], colorByCols = colorByCols)
            ggplot2::ggsave(filename=paste0(plotDir, "annoPCPlot","_", PCSTOANNOTATE[1],"_", PCSTOANNOTATE[3], "_", inputID, ".pdf"),
                            plot = colorPlot, device = "pdf",
                            limitsize=FALSE)
            colorPlot = colorClusterPlots(clusteredDF = cbind(pcScores[sharedSamples, ], pMeta[sharedSamples, ]),
                                          plotCols = PCSTOANNOTATE[c(2,3)], colorByCols = colorByCols)
            ggplot2::ggsave(filename=paste0(plotDir, "annoPCPlot","_", PCSTOANNOTATE[2],"_", PCSTOANNOTATE[3], "_", inputID, ".pdf"),
                            plot = colorPlot, device = "pdf",
                            limitsize=FALSE)
            
        }
        },error=function(cond) {message("Error with colored target variable plots")}, finally = {})
    }
    
    #########################################################################
    if (makePairwisePCPlots) {
        # # make plots of pairwise PC scores
        # subPCA = dimRedOnRS(regionSet = GRList[[regionSetName[i]]], methylData = bigSharedC$methylProp, mCoord = bigSharedC$coordinates, drMethod = "pca")
        # pcaWithAnno = cbind(subPCA$x, patientMetadata_pqc)
        # plotPairwiseColPCs(pcsToPlot = 1:7, pcaWithAnno = pcaWithAnno, 
        #                    colorByCols = colorByCols["NPM1"], getPlotSubdir(paste0(regionSetID[i], "/")), nameString = regionSetID[i])
        # 
    }
    
}
