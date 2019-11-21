
# for rsEnrichment

###############################################################################
# loading data, assigning data to generic names for below script (so below 
# script can be kept the same regardless of the data)

# rsEnSortedInd
# GRList
# coordinateDT
# loadingMat
# rsEnrichment
# mPCA
# plotSubdir # directory for this analysis. A sub directory in plotSubdir will be added for the specific data input (inputID)
# methylData
# inputID # since a cache is saved, this marks what the methylation data source was
# Sys.getenv("PLOTS") should be set

# CAPS variables should not be changed later in the script 
# (but R does not have immutable objects)
# PCSTOANNOTATE 

#### plots that will be created and script specific parameters for them 
## "comparePCHeatmap"
# PCsToAnnotate_cPCH
## "methylAlongPC"
# topRSToPlotNum
# PCsToAnnotate_mAPC = PCSTOANNOTATE[1:5]
## "regionQuantileByPC"
# PCsToAnnotate_rQBPC = c("PC1m4", "PC1p3", PCSTOANNOTATE)
# topRSInd_rQBPC = unique(unlist(rsEnSortedInd[1:10, ])) # get top 10 region sets from each PC
## "pcFromSubset Correlation Heatmap"
# PCsToAnnotate_pcFSCH
# topRSInd_pcFSCH = unique(unlist(rsEnSortedInd[1:10, ])) # get top region sets from each PC
## "region set Overlapping Cytosine Proportion" (rsOLCP)
## proportion of cytosines from region set that are shared with other region set
# topRSInd_rsOLCP = unique(unlist(rsEnSortedInd[1:10, ]))
## "meta region loading profiles" (mrLP)
# topRSInd_mrLP
# PCsToAnnotate_mrLP = PCSTOANNOTATE

# by default make all plots unless already specified by setting to FALSE
if (!exists("makeCPCH")) {
    makeCPCH = TRUE
}
if (!exists("makeMAPC")) {
    makeMAPC = TRUE
}
if (!exists("makeRQBPC")) {
    makeRQBPC = TRUE
}
if (!exists("makePCFSCH")) {
    makePCFSCH = TRUE
}
if (!exists("makerRSOLCP")) {
    makeRSOLCP = TRUE
}
if (!exists("makeMRLP")) {
    makeMRLP = TRUE
}


#################################################################################
# place to save plots
plotSubdir = paste0(plotSubdir, "/")
if (!dir.exists(paste0(Sys.getenv("PLOTS"), plotSubdir))) {
    dir.create(paste0(Sys.getenv("PLOTS"), plotSubdir), recursive = TRUE)
}

inputID = addUnderscore(inputID, side = "left")

thisPlotSubdir = paste0(plotSubdir, "COCOA_plots", inputID, "/")
createPlotSubdir(thisPlotSubdir)



source(paste0(Sys.getenv("CODE"), "aml_e3999/src/00-genericFunctions.R"))
source(paste0(Sys.getenv("CODE"), "PCRSA_extra/R/visualization.R"))
source(paste0(Sys.getenv("CODE"), "PCRSA_extra/R/analysis.R"))




library(grid)
library(ggplot2)
library(COCOA)
library(ComplexHeatmap)
devtools::load_all(ffCode("COCOA/"))

##################################################################################
# comparePCHeatmap
# visualization of enrichment score results across PCs
# see if region set high in one PC is also high in others


# number of plots = length(PCsToRankBy). one plot for each
# multiple plots in one pdf
# TODO: filter out low coverage region sets
# for rsEnrichment
if (makeCPCH) {
    comparePCHeatmap(rsScores=rsEnrichment, 
                     PCsToRankBy=PCsToAnnotate_cPCH, 
                     PCsToInclude=PCsToAnnotate_cPCH,
                     fileName=paste0(Sys.getenv("PLOTS"), thisPlotSubdir, 
                                     "rsEnrichHeatmap", inputID, ".pdf"))
}



##################################################################################
# methylAlongPC
# looking at methylation level data at individual cytosines ordered by PC 
# only looking at regions with high average loading scores
# still individual cytosine methylation

if (makeMAPC) {
    # one pdf for each PC given.
    for (i in seq_along(PCsToAnnotate_mAPC)) {
        
        # top region sets for this PC
        rsInd = as.numeric(as.matrix(rsEnSortedInd[1:topRSToPlotNum, PCsToAnnotate_mAPC[i]])) # original index
        
        grDevices::pdf(paste0(Sys.getenv("PLOTS"), thisPlotSubdir, "regionMethylHeatmaps", PCsToAnnotate_mAPC[i], inputID, ".pdf"), width = 11, height = 8.5 * topRSToPlotNum)
        
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
}

###################################################################################
# regionQuantileByPC
# comparing loading scores/percentiles for individual regions among PCs
# need region sets and PCA loadings

if (makeRQBPC) {
    if (!dir.exists(paste0(Sys.getenv("PLOTS"), thisPlotSubdir, "regionByPC/"))) {
        dir.create(paste0(Sys.getenv("PLOTS"), thisPlotSubdir, "regionByPC/"))
    }
    
    grDevices::pdf(paste0(Sys.getenv("PLOTS"), thisPlotSubdir, 
                          "regionByPC/regionPercentileByPC", inputID, ".pdf"), 
                   width = 11, height = 8.5 * length(topRSInd_rQBPC))
    
    ## ranking in terms of percentiles in case there were different distributions of loading scores for each PC
    # if there are too many regions, will try to cluster and cause memory error:
    # cannot allocate vector of size X Gb,
    # fix this by decreasing maxRegionsToPlot or use cluster_rows=FALSE
        multiRegionQuantileByPC(signal=loadingMat, signalCoord=coordinateDT, 
                           GRList=GRList[topRSInd_rQBPC], 
                           rsNames=paste0(rsEnrichment$rsName[topRSInd_rQBPC], " : ", rsEnrichment$rsDescription[topRSInd_rQBPC]), 
                            signalCol=PCsToAnnotate_rQBPC, maxRegionsToPlot = 5000,
                           cluster_rows = TRUE, absVal = TRUE)
    
    dev.off()
}

##################################################################################
# pcFromSubset Heatmap
# seeing whether a subset of cytosines (ie a single region set) can
# recapitulate PC score from all cytosines

if (makePCFSCH) {
    .regionSetList = GRList[topRSInd_pcFSCH] 
    regionSetName = paste0(rsEnrichment$rsName[topRSInd_pcFSCH], " : ", rsEnrichment$rsDescription[topRSInd_pcFSCH])
    names(.regionSetList) <-  regionSetName
    subsetCorList = lapply(X = as.list(.regionSetList), FUN = function(x) pcFromSubset(regionSet = x, 
                                                                                       mPCA = mPCA, 
                                                                                       methylData = methylData, 
                                                                                       mCoord = coordinateDT, 
                                                                                       pc = PCsToAnnotate_pcFSCH,
                                                                                       returnCor = TRUE))
    subsetCorMat = do.call(rbind, subsetCorList)
    colnames(subsetCorMat) <- PCsToAnnotate_pcFSCH
    
    simpleCache(paste0("subsetCorMat", inputID), {subsetCorMat}, recreate=TRUE)
    
    grDevices::pdf(paste0(Sys.getenv("PLOTS"), thisPlotSubdir, "subsetCorRSbyPC", inputID, ".pdf"), width = 8.5, 11)
    
    draw(Heatmap(matrix = subsetCorMat, col = c("black", "orange"), cluster_rows = FALSE, cluster_columns = FALSE, 
            column_title = , cell_fun = function(.j, .i, x, y, width, height, fill, mat=subsetCorMat) {
                grid.text(sprintf("%.2f", mat[.i, .j]), x, y, gp = gpar(fontsize = 10))
            }))
    dev.off()
    
    
    # plotting correlation between a PC and the "PC-subset" score 
    # derived from only loading values of CpGs within a certain region set 
    # do this for top few region sets (i in the loop) for each PC
    for (j in seq_along(PCsToAnnotate_pcFSCH)) {
        grDevices::pdf(paste0(Sys.getenv("PLOTS"), 
                              thisPlotSubdir, "/subsetPC_PC_Scatter_", 
                              PCsToAnnotate_pcFSCH[j], addUnderscore(inputID, side="left"), ".pdf"))
        for (i in 1:10) {
            # returns a "PC-subset" score for each sample
            subScores = pcFromSubset(regionSet = GRList[as.numeric(as.data.frame(rsEnSortedInd)[i, PCsToAnnotate_pcFSCH[j]])], 
                                     mPCA = mPCA, 
                                     methylData = methylData, 
                                     mCoord = coordinateDT, 
                                     pc = PCsToAnnotate_pcFSCH[j],
                                     returnCor = FALSE)    
            
            plot(x = subScores, y = mPCA$x[, PCsToAnnotate_pcFSCH[j]], xlab= "Scores from only CpGs in this region set",
                 ylab = PCsToAnnotate_pcFSCH[j], main = paste0(as.data.frame(rsEnrichment)$rsName[as.numeric(rsEnSortedInd[i, PCsToAnnotate_pcFSCH[j]])], 
                                                       " : ", as.data.table(rsEnrichment)$rsDescription[as.numeric(rsEnSortedInd[i, PCsToAnnotate_pcFSCH[j]])]))   
            
        }
        dev.off()
    }
}


################################################################################
# seeing how much overlap there is between region sets
# based on overlap of covered cytosines, not the regions themselves

if (makeRSOLCP) {
    # total regions in column region sets are the denominator for the proportion
    .regionSetList = GRList[topRSInd_rsOLCP] 
    pOL = percentCOverlap(mCoord = MIRA:::dtToGr(coordinateDT), 
                          GRList = .regionSetList) 
    
    grDevices::pdf(paste0(Sys.getenv("PLOTS"), thisPlotSubdir, "topRSOverlap", inputID, ".pdf"), width = 25, height = 25)
    
    draw(Heatmap(matrix = pOL[[1]], col = c("black", "yellow"), cluster_rows = FALSE, cluster_columns = FALSE, 
            cell_fun = function(j, i, x, y, width, height, fill, mat=pOL[[1]]) {
                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
            }))
    
    # numbers not included on plot squares
    # Heatmap(matrix = pOL[[1]], col = c("gray14", "gold"), cluster_rows = FALSE, cluster_columns = FALSE)
    dev.off()
}

####################################################################
# "meta region loading profiles" (mrLP)
# meta-region plots of region surrounding regions in region set
# check whether enrichment is specific to this region set by
# seeing if loading values have a spike in the center of these region sets
# compared to surrounding genome 

if (makeMRLP) {
    .regionSetList = GRList[topRSInd_mrLP]
    .regionSetList = lapply(.regionSetList, resize, width = 14000, fix="center")
    .rsNames = paste0(rsEnrichment$rsName[topRSInd_mrLP], " : ", rsEnrichment$rsDescription[topRSInd_mrLP])
    
    # returns a list: one item is grob, one item is list of binned data tables
    mrPlotOutput = makeMetaRegionPlots(signal = loadingMat, signalCoord = coordinateDT, GRList = .regionSetList, 
                        rsNames = .rsNames, signalCol = PCsToAnnotate_mrLP, binNum = 21, 
                        aggrMethod = "default")
    pcProf = mrPlotOutput$metaRegionData

    simpleCache(paste0("pcProf14k", inputID), {
        pcProf
    }, recreate = TRUE)
   
    ggsave(filename = paste0(Sys.getenv("PLOTS"), thisPlotSubdir,
                             "/metaRegionLoadingProfiles", 
                             inputID, ".pdf"), plot = mrPlotOutput$grob, device = "pdf", limitsize = FALSE)
    
    # check PPARG.bed, Jaspar motifs (had 18 rows instead of 21)
}

#####################################################################


