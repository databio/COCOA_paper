# load libraries and prepare environment for other scripts

library(dplyr)
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
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-dataProcessingFunctions.R" ))
library(MultiAssayExperiment)
library(folderfun)

# source(paste0(Sys.getenv("CODE"), "COCOA/R/COCOA.R"))
# source(paste0(Sys.getenv("CODE"), "COCOA/R/visualization.R"))



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
devtools::load_all(ffCode("COCOA/"))

createPlotSubdir <- function(plotSubdir) {
    if (!dir.exists(ffPlot(plotSubdir))) {
        dir.create(ffPlot(plotSubdir))
    }
}

# for ggplot2
theme_set(theme_classic() + 
              theme(axis.text = element_text(colour = "black", size = 15), axis.ticks = element_line(colour = "black")))

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

# @param dataMat columns of dataMat should be samples/patients, rows should be genomic signal
# (each row corresponds to one genomic coordinate/range)
# @param featureMat Rows should be samples, columns should be "features" 
# (whatever you want to get correlation with: eg PC scores),
# all columns in featureMat will be used (subset when passing to function
# in order to not use all columns)
# @param centerDataMat logical object. Should rows in dataMat be centered based on
# their means? (subtracting row means from each row)
# @param centerFeatureMat. logical.
# @param testType character object. Can be "cor" (Pearson correlation),
# "spearmanCor (Spearman correlation)
# "pcor" (partial correlation), "cov" (covariance (Pearson)), "spearmanCov" 
# (Spearman covariance)
# @param covariate
#
# If a row in dataMat has 0 stand. deviation, correlation will be set to 0
# instead of NA as would be done by cor()
#
# returns a matrix where rows are the genomic signal (eg a CpG or region) and
# columns are the columns of featureMat
# @examples dataMat = matrix(rnorm(50), 5, 10)
# featureMat = matrix(rnorm(20), 10, 2)
createCorFeatureMat = function(dataMat, featureMat, 
                               centerDataMat=TRUE, centerFeatureMat = TRUE, 
                               testType="cor", covariate=NULL) {
    
    featureMat = as.matrix(featureMat)
    featureNames = colnames(featureMat)
    nFeatures = ncol(featureMat)
    nDataDims = nrow(dataMat)
    
    if (centerDataMat) {
        cpgMeans = rowMeans(dataMat, na.rm = TRUE)
        # centering before calculating correlation
        dataMat = apply(X = dataMat, MARGIN = 2, function(x) x - cpgMeans)
        
    }
    
    if (centerFeatureMat) {
        featureMeans = colMeans(featureMat, na.rm = TRUE)
        # centering before calculating correlation
        featureMat = t(apply(X = t(featureMat), MARGIN = 2, function(x) x - featureMeans))
        if (dim(featureMat)[1] == 1) {
            featureMat = t(featureMat)
        }
        featureMat = as.matrix(featureMat)
    }
    
    dataMat = data.table::copy(as.data.frame(t(dataMat)))
    
    
    if (testType == "cor") {
        # create feature correlation matrix with PCs (rows: features/CpGs, columns:PCs)
        # how much do features correlate with each PC?
        
        # featurePCCor = as.data.frame(matrix(rep(0, nFeatures * nDataDims), nrow=nDataDims, ncol=nFeatures))
        # for (i in 1:nFeatures) {
        #     for (j in 1:nDataDims) {
        #         featurePCCor[j, i] = cor(x = featureMat[, i], y = dataMat[, j], use="pairwise.complete.obs")
        #     }
        #     
        # }
        featurePCCor = apply(X = featureMat, MARGIN = 2, function(y) apply(X = dataMat, 2,
                                                                           FUN = function(x) cor(x = x, y,
                                                                                                 use="pairwise.complete.obs",
                                                                                                 method="pearson")))
    } else if (testType == "spearmanCor") {
        featurePCCor = apply(X = featureMat, MARGIN = 2, function(y) apply(X = dataMat, 2,
                                                                           FUN = function(x) cor(x = x, y,
                                                                                                 use="pairwise.complete.obs",
                                                                                                 method="spearman")))
    } else if (testType == "pcor") {
        featurePCCor = apply(X = featureMat, MARGIN = 2, function(y) apply(X = dataMat, 2, 
                                                                           FUN = function(x) pcor.test(x = x, y=y,
                                                                                                       z=covariate)$estimate))
        
    } else if (testType == "cov") {
        # featurePCCor = as.data.frame(matrix(rep(0, nFeatures * nDataDims), nrow=nDataDims, ncol=nFeatures))
        # for (i in 1:nFeatures) {
        #     for (j in 1:nDataDims) {
        #         featurePCCor[j, i] = cov(x = featureMat[, i], y = dataMat[, j], use="pairwise.complete.obs")
        #     }
        #     
        # }
        
        featurePCCor = apply(X = featureMat, MARGIN = 2, function(y) apply(X = dataMat, 2,
                                                                           FUN = function(x) cov(x = x, y,
                                                                                                 use="pairwise.complete.obs")))
    } else {
        stop("invalid testType")
    }
    
    
    # if standard deviation of the data was zero, NA will be produced
    # set to 0 because no standard deviation means no correlation with attribute of interest
    featurePCCor[is.na(featurePCCor)] = 0
    colnames(featurePCCor) <- featureNames
    
    return(featurePCCor)
    # corLoadRatio = signal[, signalCol] / featurePCCor 
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

#' @param nullDistDF a data.frame. Has null distributions for a single 
#' region set. Each column corresponds to a null distribution for that 
#' region set for a given variable/sample attribute.   
#' 
#' @return Returns a list. Each list item is a "fitdist" object which is 
#' a fitted function 
#' for one of the columns in nullDistDF (output of fitdist() from
#' fitdistrplus.
#' These are gamma distributions and can be used to get p values for the 
#' null distribution so that a large number of permutations 
#' are not required.The list is in order of the columns and will
#' have the names of the data.frame columns. 

fitGammaNullDistr <- function(nullDistDF, method="mme", force=FALSE) {
    if (force) {
        modelList = apply(X = nullDistDF, 
                          MARGIN = 2, 
                          FUN = function(x) fitdistrplus::fitdist(data=x, 
                                                                            distr = "gamma", 
                                                                            method=method))
    } else {
        # try "method". if it fails, do method = "mme"
        modelList = apply(X = nullDistDF, 
                          MARGIN = 2, 
                          FUN = function(x) tryCatch({fitdistrplus::fitdist(data=x, 
                                                                            distr = "gamma", 
                                                                            method=method)}, 
                                                     error = function(e) {fitdistrplus::fitdist(data=x, 
                                                                                                 distr = "gamma", 
                                                                                                 method="mme")}))   
    }
    
    return(modelList)    
}

#' Get p value after fitting a gamma distribution to the null distribution
#' @param nullDistList list of a data.frame. Each list item 
#' has null distributions for a single 
#' region set. Each column corresponds to a null distribution for that 
#' region set for a given variable/sample attribute.   
#' @param scores a data.frame. Has same columns as nullDistDF. One row per
#' region set (should be in same order as nullDistDF) The scores
#' that will be used to get p values.
#' @param realScoreInDist logical. Should the actual score (from 
#' test with no permutations) be included in the null distribution 
#' when fitting the gamma distribution
#' 
#' @return Returns a data.frame with p values, one column for each col in
#' scores and nullDistDF 

getGammaPVal <- function(scores, nullDistList, method="mme", realScoreInDist=FALSE, force=FALSE) {
    
    # make sure the same columns are present/in the same order
    
    
    colsToAnnotate = colnames(scores)
    
    if (realScoreInDist) {
        # to get a more accurate gamma distribution, include the score from unpermuted test.
        # add to each null distribution
        for (i in 1:nrow(scores)) {
            nullDistList[[i]] = rbind(nullDistList[[i]], as.numeric(scores[i, ]))
        }
    }
 
    
    # returns list, each item in list is also a list.
    # in sub list: each col in nullDistDF has one list item
    fittedDistList = lapply(X = nullDistList, function(x) fitGammaNullDistr(nullDistDF = x[, colsToAnnotate, drop=FALSE], 
                                                                            method=method, 
                                                                            force=force))
    
    pValList = list()
    # once for each region set
    for (i in seq_along(nullDistList)) {
        
        pValList[[i]] = as.data.frame(t(pGammaList(scoreVec = as.numeric(scores[i, ]), 
                                                   fitDistrList = fittedDistList[[i]])))
    }
    pValDF = as.data.frame(rbindlist(pValList))
    
    
    return(pValDF)
}

pGammaList <- function(scoreVec, fitDistrList) {
    pValVec = mapply(FUN = function(x, y) pgamma(q = y, 
                      shape = x$estimate["shape"], 
                      rate = x$estimate["rate"], 
                      lower.tail = FALSE), x = fitDistrList, y = scoreVec)
    return(pValVec)
}

#####################################################################
# functions for parallel permutations

# @param genomicSignal columns of dataMat should be samples/patients, rows should be genomic signal
# (each row corresponds to one genomic coordinate/range)
# @param sampleLabels Matrix or data.frame. Rows should be samples, 
# columns should be "features" 
# (whatever you want to get correlation with: eg PC scores),
# all columns in featureMat will be used (subset when passing to function
# in order to not use all columns)
# @param calcCols character. the columns in `sampleLabels` for which to calculate
# correlation and then to run COCOA on
corPerm <- function(randomInd, genomicSignal, 
                    signalCoord, GRList, calcCols,
                    sampleLabels, variationMetric = "cor") {
    
    # if vector is given, return error
    if (is.null(dim(sampleLabels))) {
        stop("`sampleLabels` should be a matrix or data.frame")
    }
    
    if (any(!(calcCols %in% colnames(sampleLabels)))) {
        stop("Not all specified columns are present in `sampleLabels`")
    }
    
    # subset to only calcCols
    sampleLabels = sampleLabels[, calcCols, drop=FALSE]
    
    # because names are dropped for a single column data.frame when indexing
    # single col data.frame is automatically converted to numeric
    featureNames = colnames(sampleLabels)
    # reorder the sample labels
    sampleLabels = data.frame(sampleLabels[randomInd, ])
    colnames(sampleLabels) = featureNames
    
    # calculate correlation
    featureLabelCor = createCorFeatureMat(dataMat = genomicSignal, 
                                          featureMat = sampleLabels, 
                                          centerDataMat = TRUE, 
                                          centerFeatureMat = TRUE,
                                          testType = variationMetric)
    
    # run COCOA
    thisPermRes = runCOCOA(signal=featureLabelCor, 
                           signalCoord=signalCoord, GRList=GRList, 
                           signalCol = calcCols, 
                           scoringMetric = "default", verbose = TRUE)
    
    # return
    return(thisPermRes)
    
}
# This function will take a list of results of permutation tests that included
# many region sets and return a data.frame/data.table with the null
# distribution for a single region set (row)
# @param resultsList each item in the list is a data.frame, one item for
# each permutation with the results of that permutation. Each row in the 
# data.frame is a region set. Rows in all the data.frames should be
# in the same order.
# @param rsInd numeric. The row number for the region set of interest.
extractNullDist <- function(resultsList, rsInd) {
    rowList = lapply(resultsList, FUN = function(x) x[rsInd, ])
    rsNullDist = as.data.frame(rbindlist(rowList))
    return(rsNullDist)
}


################## not sure whether to add this one to COCOA
# convenience function to quickly make meta-region loading profiles
# first calculates profiles, then normalizes and plots
# output is marrangeGrob and should be saved with ggsave
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

makeMetaRegionPlots <- function(signal, signalCoord, GRList, rsNames, signalCol, binNum, aggrMethod="default") {
    
    pcProf = lapply(X = GRList, function(x) getMetaRegionProfile(signal = signal, 
                                                              signalCoord = signalCoord, 
                                                              regionSet = x, signalCol = signalCol,
                                                              binNum = binNum, aggrMethod=aggrMethod))
    
    pcP = copy(pcProf)
    # this check must be done before converting items of list to data.table
    # otherwise NULL will be converted to a data.table
    notNull = !vapply(X = pcP, FUN = is.null, FUN.VALUE = TRUE)
    pcP = pcP[notNull]
    rsNames = rsNames[notNull]
    pcP = lapply(pcP,FUN = as.data.table)
    
    
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
    
    # average loading value from each PC to normalize so PCs can be compared with each other
    avLoad = apply(X = signal[, signalCol], MARGIN = 2, FUN = function(x) mean(abs(x)))
    
    # normalize
    # pcP = lapply(pcP, FUN = function(x) t(apply(X = x, MARGIN = 1, FUN = function(y) y - c(0, avLoad))))
    # for (i in seq_along(signalCol)) {
    #     # by reference
    #     lapply(pcP, FUN = function(x) x[, 1])
    # }
    pcP = lapply(pcP, FUN = function(x) x[, mapply(FUN = function(y, z) get(y) - z, y=signalCol, z = avLoad)])
    pcP = lapply(pcP, FUN = function(x) data.table(regionGroupID=1:nrow(x), x))
    
    # for the plot scale
    maxVal = max(sapply(pcP, FUN = function(x) max(x[, .SD, .SDcols=signalCol])))
    minVal = min(sapply(pcP, FUN = function(x) min(x[, .SD, .SDcols=signalCol])))
    
    # convert to long format for plots
    pcP = lapply(X = pcP, FUN = function(x) tidyr::gather(data = x, key = "PC", value="loading_value", signalCol))
    pcP = lapply(X = pcP, as.data.table)
    pcP = lapply(pcP, function(x) x[, PC := factor(PC, levels = signalCol)])
    
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
    
    metaRegionProfileInfo = list(multiProfileP, pcProf)
    setattr(metaRegionProfileInfo, "names", c("grob", "metaRegionData"))
    
    return(metaRegionProfileInfo)
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

# from dev version of COCOA
signalAlongPC <- function(genomicSignal, signalCoord, regionSet, 
                          sampleScores, orderByCol="PC1", topXVariables=NULL, 
                          variableScores=NULL,
                          cluster_columns = FALSE, 
                          cluster_rows = FALSE, 
                          row_title = "Sample",
                          column_title = "Genomic Signal", 
                          column_title_side = "bottom",
                          name = "Genomic Signal Value",
                          col = c("blue", "#EEEEEE", "red"), ...) {
    
    
    if (!(is(genomicSignal, "matrix") || is(genomicSignal, "data.frame"))) {
        stop("genomicSignal should be a matrix or data.frame. Check object class.")
    }
    
    # test for appropriateness of inputs/right format
    if (is(signalCoord, "GRanges")) {
        coordGR <- signalCoord
    } else if (is(signalCoord, "data.frame")) {
        # UPDATE: does the work on data.frames that are not data.tables?
        coordGR <- dtToGr(signalCoord)
    } else {
        stop("signalCoord should be a data.frame or GRanges object.")
    }
    
    if (!(is(sampleScores, "matrix") || is(sampleScores, "data.frame"))) {
        stop("sampleScores should be a matrix or data.frame.")
    }
    
    
    # PCA object must have subject_ID as row.names (corresponding 
    # to column names of genomicSignal)
    if (sum(row.names(sampleScores) %in% colnames(genomicSignal)) < 2) {
        stop(cleanws("Sample names on pca data (row names) 
                     must match sample names on methylation
                     (column names)"))
    }
    
    
    if (!is(regionSet, "GRanges")) {
        stop("regionSet should be a GRanges object. Check object class.")
    }
    
    if (!is.null(topXVariables)) {
        if (is.null(variableScores)) {
            stop("To plot the topXVariables, variableScores must be given.")
        }
        if (!(length(variableScores) == nrow(genomicSignal))) {
            stop("length(variableScores) should equal nrow(genomicSignal)")
        }
        
    }
    
    
    
    # coordGR =
    olList <- findOverlaps(regionSet, coordGR)
    # regionHitInd <- sort(unique(queryHits(olList)))
    cytosineHitInd <- sort(unique(subjectHits(olList)))
    thisRSMData <- t(genomicSignal[cytosineHitInd, ])
    nRegion = length(unique(queryHits(olList)))
    # get top variables
    if (!is.null(topXVariables)) {
        if (nRegion > topXVariables) {
            # select top variables
            thisRSMData <- thisRSMData[, order(variableScores[cytosineHitInd], decreasing = TRUE)][, 1:topXVariables]
        }
    }
    subject_ID <- row.names(thisRSMData)
    # centeredPCAMeth <- t(apply(t(genomicSignal), 1, 
    #                            function(x) x - pcaData$center)) #center first 
    # reducedValsPCA <- centeredPCAMeth %*% pcaData$rotation
    # reducedValsPCA <- pcaData$x
    # pcaData must have subject_ID as row name
    thisRSMData <- thisRSMData[names(sort(sampleScores[, orderByCol], 
                                          decreasing = TRUE)), ]
    
    message(paste0("Number of cytosines: ", length(cytosineHitInd)))
    message(paste0("Number of regions: ", nRegion))
    

    
    ComplexHeatmap::Heatmap(thisRSMData, 
                            col = col,
                            row_title = row_title,
                            column_title = column_title,
                            column_title_side = column_title_side,
                            cluster_rows = cluster_rows, 
                            cluster_columns = cluster_columns, 
                            name = name, ...)
}

