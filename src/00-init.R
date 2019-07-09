# load libraries and prepare environment for other scripts

library(LOLA)
library(simpleCache)
library(data.table)
library(ggplot2)
library(GenomicRanges) # GRangesList, resize
# library(caret)
# library(RGenomeUtils)
library(MIRA)
library(ComplexHeatmap)
library(gridExtra) #marrangeGrob for colorClusterPlots()
# some of the environmental variables from aml/.../00-init.R will need to be reset
source(paste0(Sys.getenv("CODE"), "aml_e3999/src/00-genericFunctions.R" )) 
library(MultiAssayExperiment)
library(folderfun)

# source(paste0(Sys.getenv("CODE"), "COCOA/R/COCOA.R"))
# source(paste0(Sys.getenv("CODE"), "COCOA/R/visualization.R"))
library(COCOA)

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

process_brca_expr = function(exprDF) {
    
}


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
setff("ProjCode", paste0(Sys.getenv("CODE"), "COCOA_paper/"))
setff("Data", Sys.getenv("DATA"))
dirCode = function(.file="") {
    return(paste0(Sys.getenv("CODE"), "COCOA_paper/", .file))
}

# for ggplot2
theme_set(theme_classic() + 
              theme(axis.text = element_text(colour = "black", size = 15), axis.ticks = element_line(colour = "black")))

# set environment
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))

############# functions to load data easily ############################
loadBRCADNAm <- function(signalMat=TRUE, signalCoord=TRUE, 
                         loadingMat=TRUE, pcScores=TRUE,
                         .env=parent.frame(n=1)) {
    if (signalMat) {
        simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList")
        #restrict patients included in this analysis
        patientMetadata = patientMetadata[patientMetadata$subject_ID %in% 
                                              colnames(brcaMList[["methylProp"]]), ]
        # patientMetadata should have already screened out patients without ER/PGR status
        # resulting in 657 patients
        hasER_PGR_IDs = patientMetadata[, subject_ID]
        filteredMData = brcaMList[["methylProp"]][, hasER_PGR_IDs] 
        
        assign("signalMat", filteredMData, envir=.env)
    }
    
    if (signalCoord) {
        simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList",
                    reload = FALSE, recreate = FALSE)
        assign("signalCoord", brcaMList$coordinates, envir=.env)
    }
    
    if (loadingMat) {
        simpleCache("allMPCA_657", assignToVariable = "allMPCA")
        loadingMat = allMPCA$rotation
        assign("loadingMat", loadingMat, envir=.env)
    }
    if (pcScores) {
        simpleCache("allMPCA_657", assignToVariable = "allMPCA", reload = FALSE)
        pcScores = allMPCA$x
        assign("pcScores", pcScores, envir=.env)
    }
}

loadMOFAData <- function(methylMat=TRUE, signalCoord=TRUE, latentFactors=TRUE, 
                         factorWeights=FALSE, .env=parent.frame(n=1)) {
    
    library(ExperimentHub)
    library("SummarizedExperiment")
    # library(MOFAtools)
    library(MOFA)
    library("MOFAdata")
    library("FDb.InfiniumMethylation.hg19")
    library(MultiAssayExperiment)
    
    if (signalCoord) {
        cllMethyl = prepareCLLMethyl()
        methCoord = cllMethyl$methylCoord
        assign(x = "signalCoord", value = methCoord, envir = .env)
        
    }
    
    # if both are loaded, put samples in same order
    if (methylMat || latentFactors) {
        
        # load but don't yet assign
        if (methylMat) {
            cllMethyl = prepareCLLMethyl()
            methData = cllMethyl$methylProp
        }
        
        # load but don't yet assign
        if (latentFactors) {
            # Loading an existing trained model
            filepath <- system.file("extdata", "CLL_model.hdf5",
                                    package = "MOFAdata")
            
            MOFAobject <- loadModel(filepath)
            
            latentFactors <- getFactors(
                MOFAobject,
                as.data.frame = FALSE
            )
        }
        if (methylMat && latentFactors) {
            sharedNames = colnames(methData)[colnames(methData) %in% 
                                                 row.names(latentFactors)]
            latentFactors = latentFactors[sharedNames, ]
            methData = methData[, sharedNames]
            
            assign("latentFactors", latentFactors, envir = .env)
            assign("methylMat", methData, envir = .env)
        } else if (methylMat) {
            assign("methylMat", methData, envir = .env)
        } else if (latentFactors) {
            assign("latentFactors", latentFactors, envir = .env)
        }
    }
    
    if (factorWeights) {
        # Loading an existing trained model
        filepath <- system.file("extdata", "CLL_model.hdf5",
                                package = "MOFAdata")
        
        MOFAobject <- loadModel(filepath)
        
        MOFAweights <- getWeights(
            MOFAobject, 
            views = "Methylation", 
            factors = "all", 
            as.data.frame = TRUE
        )
        assign("factorWeights", MOFAweights, envir = .env)
    }
    
}

# loads GRList, rsName, rsDescription

loadGRList <- function(genomeV = "hg38", .env=parent.frame(n=1)) {
    if (genomeV == "hg38") {
        source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))
        names(GRList) = rsName
        assign(x = "GRList", GRList, envir = .env)
        assign(x = "rsName", rsName, envir = .env)
        assign(x = "rsDescription", rsDescription, envir = .env)
    } else if (genomeV == "hg19") {
        source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regionDB_hg19.R"))
        names(GRList) = rsName
        assign(x = "GRList", GRList, envir = .env)
        assign(x = "rsName", rsName, envir = .env)
        assign(x = "rsDescription", rsDescription, envir = .env)
    } else {
        stop("Only hg38 and hg19 are available in this function currently.")
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
    rsRankInd = rsRankingIndex(rsScores=rsScores, PCsToAnnotate=scoreColName)
    
    
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


################################################################################

# @param dataMat columns of dataMat should be samples/patients, rows should be genomic signal
# (each row corresponds to one genomic coordinate/range)
# @param featureMat Rows should be samples, columns should be "features" 
# (whatever you want to get correlation with: eg PC scores),
# all columns in featureMat will be used (subset when passing to function
# in order to not use all columns)
# @param center logical object. Should rows in dataMat be centered based on
# their means? (subtracting row means from each row)
#
# If a row in dataMat has 0 stand. deviation, correlation will be set to 0
# instead of NA as would be done by cor()
#
# returns a matrix where rows are the genomic signal (eg a CpG or region) and
# columns are the columns of featureMat
createCorFeatureMat = function(dataMat, featureMat, 
                               centerDataMat=TRUE, centerFeatureMat = TRUE) {
    
    featureMat = as.matrix(featureMat)
    
    if (centerDataMat) {
        cpgMeans = rowMeans(dataMat)
        # centering before calculating correlation
        dataMat = apply(X = dataMat, MARGIN = 2, function(x) x - cpgMeans)
        
    }
    
    if (centerFeatureMat) {
        featureMeans = colMeans(featureMat)
        # centering before calculating correlation
        featureMat = t(apply(X = featureMat, MARGIN = 1, function(x) x - featureMeans))
        if (dim(featureMat)[1] == 1) {
            featureMat = t(featureMat)
        }
        featureMat = as.matrix(featureMat)
    }
    
    
    # create feature correlation matrix with PCs (rows: features/CpGs, columns:PCs)
    # how much do features correlate with each PC?
    featurePCCor = apply(X = featureMat, MARGIN = 2, function(y) apply(X = dataMat, 1, FUN = function(x) cor(x = x, y)))
    
    # if standard deviation of the data was zero, NA will be produced
    # set to 0 because no standard deviation means no correlation with attribute of interest
    featurePCCor[is.na(featurePCCor)] = 0
    colnames(featurePCCor) <- colnames(featureMat)
    
    return(featurePCCor)
    # corLoadRatio = loadingMat[, PCsToAnnotate] / featurePCCor 
    # hist(corLoadRatio[, "PC10"])
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

# @param rsScores is a data.frame of 
# @param nullDistList list. one item per region set. Each item is a 
# data.frame with the 
# null distribution for a single region set. Each column in the data.frame
# is for a single variable (e.g. PC or latent factor)
# @param testType character. "greater", "lesser", "two-sided" Whether to
# create p values based on one sided test or not.
# @param whichMetric character. Can be "pval" or "zscore"

getPermStat <- function(rsScores, nullDistList, calcCols, testType="greater", whichMetric = "pval") {
    
    if (is(rsScores, "data.table")) {
        rsScores = as.data.frame(rsScores)
    }
    
    # get p values for a single region set (can get p val for multiple columns)
    # @param rsScore a row of values for a single region set. One 
    # value for each calcCols
    getPermStatSingle <- function(rsScore, nullDist, 
                                  calcCols, testType="greater", whichMetric = "pval") {
        
        if (is(nullDist, "data.table")) {
            nullDist = as.data.frame(nullDist)
        }
        if (whichMetric == "pval") {
            
            pVal = rep(-1, length(rsScore))
            if (testType == "greater") {
                
                for (i in seq_along(pVal)) {
                    # only for one sided test (greater than)
                    pVal[i] = 1 - ecdf(x = nullDist[, calcCols[i]])(rsScore[i])
                }
            }
            
            # (-abs(x-0.5) + 0.5) * 2
            
            thisStat = pVal
        }
        
        if (whichMetric == "zscore") {
            
            
            zScore = rep(NA, length(rsScore))
            for (i in seq_along(zScore)) {
                # only for one sided test (greater than)
                zScore[i] = (rsScore[i] - mean(nullDist[, calcCols[i]])) / sd(x = nullDist[, calcCols[i]])
            }
            
            thisStat = as.numeric(zScore)
        }
        
        
        return(thisStat)
    }
    
    
        # do once for each region set
        thisStatList = list()
        for (i in 1:nrow(rsScores)) {
            thisStatList[[i]] = as.data.frame(t(getPermStatSingle(rsScore=rsScores[i, calcCols], 
                                                              nullDist = nullDistList[[i]],
                                                              calcCols = calcCols,
                                                              whichMetric=whichMetric)))
            colnames(thisStatList[[i]]) <- calcCols
        }
        thisStat = rbindlist(thisStatList)
        # add back on annotation info
        thisStat = cbind(thisStat, rsScores[, colnames(rsScores)[!(colnames(rsScores) %in% calcCols)]])
        
        # pVals = mapply(FUN = function(x, y) getPermPvalSingle(rsScore=x, 
        #                                            nullDist = y,
        #                                            calcCols = calcCols), 
        #        x = rsScores[, calcCols], y=nullDistList)
        

    return(thisStat)
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

################## not sure whether to add this one to COCOA
# convenience function to quickly make meta-region loading profiles
# first calculates profiles, then normalizes and plots
# output is marrangeGrob and should be saved with ggsave
makeMetaRegionPlots <- function(loadingMat, signalCoord, GRList, rsNames, PCsToAnnotate, binNum, overlapMethod="single") {
    
    pcProf = lapply(X = GRList, function(x) getLoadingProfile(loadingMat = loadingMat, 
                                                              signalCoord = signalCoord, 
                                                              regionSet = x, PCsToAnnotate = PCsToAnnotate,
                                                              binNum = binNum, overlapMethod=overlapMethod))
    
    pcP = copy(pcProf)
    pcP = lapply(pcP,FUN = as.data.table)
    notNull = !vapply(X = pcP, FUN = is.null, FUN.VALUE = TRUE)
    pcP = pcP[notNull]
    rsNames = rsNames[notNull]
    
    
    # average loading value from each PC to normalize so PCs can be compared with each other
    avLoad = apply(X = loadingMat[, PCsToAnnotate], MARGIN = 2, FUN = function(x) mean(abs(x)))
    
    # normalize
    # pcP = lapply(pcP, FUN = function(x) t(apply(X = x, MARGIN = 1, FUN = function(y) y - c(0, avLoad))))
    pcP = lapply(pcP, FUN = function(x) x[, mapply(FUN = function(y, z) get(y) - z, y=PCsToAnnotate, z = avLoad)])
    pcP = lapply(pcP, FUN = function(x) data.table(regionGroupID=1:nrow(x), x))
    
    # for the plot scale
    maxVal = max(sapply(pcP, FUN = function(x) max(x[, .SD, .SDcols=PCsToAnnotate])))
    minVal = min(sapply(pcP, FUN = function(x) min(x[, .SD, .SDcols=PCsToAnnotate])))
    
    # convert to long format for plots
    pcP = lapply(X = pcP, FUN = function(x) tidyr::gather(data = x, key = "PC", value="loading_value", PCsToAnnotate))
    pcP = lapply(X = pcP, as.data.table)
    pcP = lapply(pcP, function(x) x[, PC := factor(PC, levels = PCsToAnnotate)])
    
    # stack overflow for wrapping plot title
    wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n") 
    
    
    profilePList = list()
    for (i in seq_along(pcP)) {
        
        thisRS = pcP[[i]]
        
        
        
        profilePList[[i]] = ggplot(data = thisRS, mapping = aes(x =regionGroupID , y = loading_value)) + 
            geom_line() + ylim(c(minVal, maxVal)) + facet_wrap(facets = "PC") + 
            ggtitle(label = wrapper(rsNames[i], width=30)) + xlab("Genome around Region Set, 14 kb") + 
            ylab("Normalized Loading Value") + 
            theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
        profilePList[[i]]
        # plot(pcP[[i]]$PC1, type="l") + title(rsNames[i])
        # 
        
        #xLabels = xAxisForRegionPlots2()
        
    }
    multiProfileP = marrangeGrob(profilePList, ncol = 2, nrow = 2)
    
    return(list(multiProfileP, pcProf))
}

###################### visualization ##########################

# function to make nice histograms
# add axis and overall titles
# if making multiple plots, scale/bin them consistently (consistent bins)
# able to plot multiple plots to one pdf


########## MOFA/CLL analysis #########################
# gets coordinates for CLL methyl
# returns a list with methylProp that has methylation and "methylCoord"
# that has corresponding genomic coordinates
# filters out X and Y chromosomes
prepareCLLMethyl = function(removeXY=TRUE) {
    
    # get microarray data
    eh = ExperimentHub()
    meth = eh[[names(query(eh, "CLLmethylation"))]] # EH1071
    # rows are cpgs, columns are samples
    methData = assay(meth)
    dataProbeNames = row.names(methData)
    
    
    # match probe names with genomic coordinates
    # ls('package:FDb.InfiniumMethylation.hg19')
    m450kAnno = get450k()
    length(m450kAnno)
    
    # get coordinates in same order as CLL data
    methCoord = m450kAnno[dataProbeNames]
    all(names(methCoord) == dataProbeNames)
    methCoordDT = COCOA:::grToDt(methCoord)
    # keep start coordinate as CpG site
    methCoordDT = methCoordDT[, .(chr, start)]
    
    if (removeXY) {
        xyInd = methCoordDT$chr %in% c("chrX", "chrY")
        methCoordDT = methCoordDT[!xyInd, ]
        methData = methData[!xyInd, ]
    }
    
    methCoord = COCOA:::dtToGr(methCoordDT)
    if (nrow(methCoordDT) != nrow(methData)) {
        stop("error matching probes to coordinates")
    }
    
    return(list(methylProp = methData, methylCoord = methCoord))
    
}
