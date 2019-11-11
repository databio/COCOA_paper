# load libraries and prepare environment for other scripts

# for rivanna only
.libPaths(c("/home/jtl2hk/containerFiles/3.6/", .libPaths()))

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
library(COCOA)



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
formattedCOCOAScores  <- function(scores, colsToAnnotate=paste0("PC", 1:10), numTopRS=50) {
    topRSN = numTopRS # this many top RS for each colsToAnnotate
    pRankedScores = scores
    
    pRankedScores$index = 1:nrow(pRankedScores)
    
    # get top region sets for each colsToAnnotate based on p val
    topRSZAnnoList = list()
    
    for (i in seq_along(colsToAnnotate)) {
        
        theseTopInd = dplyr::arrange(pRankedScores,
                                     desc(get(colsToAnnotate[i])))$index[1:topRSN]
        thesePValRanks = order(pRankedScores[, paste0(colsToAnnotate[i])], decreasing = FALSE)
        pRankedScores$index[thesePValRanks]
        topRSZAnnoList[[i]] = data.frame(pRankedScores[theseTopInd, c("rsName", "rsDescription", colsToAnnotate[i],
                                                                      "signalCoverage", "regionSetCoverage",
                                                                      "totalRegionNumber", "meanRegionSize")])
        
        names(topRSZAnnoList[[i]]) <- paste0(colsToAnnotate[i], "_", c("rsName", "rsDescription", "rsScore",
                                                                       "signalCoverage", "regionSetCoverage",
                                                                       "totalRegionNumber", "meanRegionSize"))
    }
    
    
    return(topRSZAnnoList)
    
    # realRSScores = rsScores
    # 
    # gPValDF = getGammaPVal(scores = realRSScores[, colsToAnnotate, drop=FALSE], nullDistList = nullDistList, method = "mme", realScoreInDist = TRUE)
    # gPValDF = apply(X = gPValDF, MARGIN = 2, FUN = function(x) p.adjust(p = x, method = correctionMethod))
    # gPValDF = cbind(gPValDF, realRSScores[, colnames(realRSScores)[!(colnames(realRSScores) %in% colsToAnnotate)]])
    # 
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
    # 
    # # get top region sets for each colsToAnnotate based on p val
    # topRSZAnnoList = list()
    # topRSN = 50 # this many top RS for each colsToAnnotate
    # for (i in seq_along(colsToAnnotate)) {
    # 
    #     theseTopInd = dplyr::arrange(pRankedScores,
    #                                  desc(get(paste0(colsToAnnotate[i], "_PValGroup"))),
    #                                  desc(get(colsToAnnotate[i])))$index[1:topRSN]
    #     thesePValRanks = order(pRankedScores[, paste0(colsToAnnotate[i], "_PVal")], decreasing = FALSE)
    #     pRankedScores$index[thesePValRanks]
    #     topRSZAnnoList[[i]] = data.frame(pRankedScores[theseTopInd, c("rsName", "rsDescription", colsToAnnotate[i],
    #                                                                   paste0(colsToAnnotate[i], "_PValGroup"),
    #                                                                   paste0(colsToAnnotate[i], "_PVal"),
    #                                                                   "signalCoverage", "regionSetCoverage",
    #                                                                   "totalRegionNumber", "meanRegionSize")])
    # 
    #     names(topRSZAnnoList[[i]]) <- paste0(colsToAnnotate[i], "_", c("rsName", "rsDescription", "rsScore",
    #                                                                    "PValGroup", "pVal", "signalCoverage", "regionSetCoverage",
    #                                                                    "totalRegionNumber", "meanRegionSize"))
    # }
    # 
    # write.csv(topRSZAnnoList, file = ffSheets(paste0("topRSPermpVals", .analysisID, ".csv")), row.names = FALSE)
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

makeMetaRegionPlots <- function(signal, signalCoord, GRList, rsNames, 
                                signalCol, binNum, 
                                # returnNormalizedVals=TRUE, 
                                aggrMethod="default", absVal=TRUE) {
    
    pcProf = lapply(X = GRList, function(x) getMetaRegionProfile(signal = signal, 
                                                              signalCoord = signalCoord, 
                                                              regionSet = x, signalCol = signalCol,
                                                              binNum = binNum, 
                                                              aggrMethod=aggrMethod,
                                                              absVal = absVal))
    
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
    if (absVal) {
        avLoad = apply(X = signal[, signalCol], MARGIN = 2, FUN = function(x) mean(abs(x)))    
    } else {
        avLoad = apply(X = signal[, signalCol], MARGIN = 2, FUN = function(x) mean(x))
    }
    
    
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
    names(profilePList) <- rsNames    

    metaRegionProfileInfo = list(multiProfileP, pcProf, profilePList)
    setattr(metaRegionProfileInfo, "names", c("grob", "metaRegionData", "plotList"))
    
    return(metaRegionProfileInfo)
}

# stack overflow for wrapping plot title
wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n") 

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
    library("cowplot")
    
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
plotAnnoScoreDist <- function(rsScores, colsToPlot, pattern, patternName=pattern) {
    
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
        rsCorP = rsCorP + geom_point(data = rsScores[as.numeric(rsScores$Group) == i, ], alpha=0.5, shape=3) 
    }
    
    return(rsCorP)
}
