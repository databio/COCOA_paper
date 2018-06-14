

# TODO make script do plots for top10 or just take that out and run script twice
# once for rsEnrichment and once for rsEnrichmentTop10

###############################################################################
# loading data, assigning data to generic names for below script (so below 
# script can be kept the same regardless of the data)
# search "script_specific" to find portions of script that may need to be 
# changed

# rsEnSortedInd
# GRList
# coordinateDT
# loadingMat
# rsEnrichment
# rsEnrichmentTop10
# mPCA
# plotSubdir # one directory per PCA. This script generates plots for one PCA.
# methylData
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


#################################################################################
# place to save plots
plotSubdir = paste0(plotSubdir, "/")
if (!dir.exists(paste0(Sys.getenv("PLOTS"), plotSubdir))) {
    dir.create(paste0(Sys.getenv("PLOTS"), plotSubdir), recursive = TRUE)
}

library(grid)
library(ComplexHeatmap)

##################################################################################
# comparePCHeatmap
# visualization of enrichment score results across PCs
# see if region set high in one PC is also high in others


# number of plots = length(PCsToRankBy). one plot for each
# TODO: filter out low coverage region sets
# for rsEnrichment
comparePCHeatmap(rsEnrichment=rsEnrichment, 
                 PCsToRankBy=PCsToAnnotate_cPCH, 
                 PCsToInclude=PCsToAnnotate_cPCH,
                 fileName=paste0(Sys.getenv("PLOTS"), plotSubdir, "rsEnrichHeatmap.pdf"))

# for rsEnrichmentTop10
comparePCHeatmap(rsEnrichment=rsEnrichmentTop10, 
                 PCsToRankBy=PCsToAnnotate_cPCH, 
                 PCsToInclude=PCsToAnnotate_cPCH,
                 fileName=paste0(Sys.getenv("PLOTS"), plotSubdir, "rsEnrichHeatmapTop10Variable.pdf"))


##################################################################################
# methylAlongPC
# looking at methylation level data at individual cytosines ordered by PC 
# only looking at regions with high average loading scores
# still individual cytosine methylation

# one pdf for each PC given.
for (i in seq_along(PCsToAnnotate_mAPC)) {
    
    # top region sets for this PC
    rsInd = as.numeric(as.matrix(rsEnSortedInd[1:topRSToPlotNum, PCsToAnnotate_mAPC[i], with=FALSE])) # original index
    
    grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir, "regionMethylHeatmaps", PCsToAnnotate_mAPC[i], ".pdf"), width = 11, height = 8.5 * topRSToPlotNum)
    
    # heatmap
    methylAlongPC(loadingMat=loadingMat, loadingThreshold=0.95,
                  pcScores=mPCA$x,
                  coordinateDT=coordinateDT,
                  methylData=methylData,
                  GRList=GRList[rsInd], orderByPC=PCsToAnnotate_mAPC[i],
                  topXRegions=50)

    # draw(Heatmap(matrix = methylData[1:1000, 1:10]))
    # plot(methylData[1:1000, 1])
    
    dev.off()
}

###################################################################################
# regionQuantileByPC
# comparing loading scores/percentiles for individual regions among PCs
# need region sets and PCA loadings


if (!dir.exists(paste0(Sys.getenv("PLOTS"), plotSubdir, "regionByPC/"))) {
    dir.create(paste0(Sys.getenv("PLOTS"), plotSubdir, "regionByPC/"))
}

grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir, 
                      "regionByPC/regionPercentileByPC", ".pdf"), 
               width = 11, height = 8.5 * length(topRSInd_rQBPC))

# ranking in terms of percentiles in case there were different distributions of loading scores for each PC
regionQuantileByPC(loadingMat=loadingMat, coordinateDT=coordinateDT, 
                   GRList=GRList[topRSInd_rQBPC], 
                   rsNames=paste0(rsEnrichment$rsName[topRSInd_rQBPC], " : ", rsEnrichment$rsDescription[topRSInd_rQBPC]), 
                   PCsToAnnotate=PCsToAnnotate_rQBPC)
dev.off()

##################################################################################
# pcFromSubset Heatmap
# seeing whether a subset of cytosines (ie a single region set) can
# recapitulate PC score from all cytosines

regionSetList = GRList[topRSInd_pcFSCH] 
regionSetName = paste0(rsEnrichment$rsName[topRSInd_pcFSCH], " : ", rsEnrichment$rsDescription[topRSInd_pcFSCH])
names(regionSetList) <-  regionSetName
subsetCorList = lapply(X = as.list(regionSetList), FUN = function(x) pcFromSubset(regionSet = x, 
                                                                                  pca = mPCA, 
                                                                                  methylData = methylData, 
                                                                                  coordinateDT = coordinateDT, 
                                                                                  PCofInterest = PCsToAnnotate_pcFSCHM,
                                                                                  returnCor = TRUE))
subsetCorMat = do.call(rbind, subsetCorList)
colnames(subsetCorMat) <- PCsToAnnotate_pcFSCHM

simpleCache("subsetCorMat", {subsetCorMat}, recreate=TRUE)

grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir, "subsetCorRSbyPC", ".pdf"), width = 8.5, 11)

# don't use i for index since it is defined as something else in cell_fun
Heatmap(matrix = subsetCorMat, cluster_rows = FALSE, cluster_columns = FALSE, 
        column_title = , cell_fun = function(j, i, x, y, width, height, fill, mat=subsetCorMat) {
            grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

################################################################################
# seeing how much overlap there is between region sets
# based on overlap of covered cytosines, not the regions themselves


# total regions in column region sets are the denominator for the proportion
regionSetList = GRList[topRSInd_rsOLCP] 
pOL = percentCOverlap(coordGR = MIRA:::dtToGr(coordinateDT), 
                      GRList = regionSetList) 

grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir, "topRSOverlap", ".pdf"), width = 25, height = 25)

Heatmap(matrix = pOL[[1]], cluster_rows = FALSE, cluster_columns = FALSE, 
        column_title = , cell_fun = function(j, i, x, y, width, height, fill, mat=pOL[[1]]) {
            grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()