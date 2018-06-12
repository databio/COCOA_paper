

###############################################################################
# loading data, assigning data to generic names for below script (so below 
# script can be kept the same regardless of the data)
# search "script_specific" to find portions of script that may need to be 
# changed

rsEnSortedInd

GRList[rsIndex] 
coordinateDT = brcaMList$coordinates
loadingMat = allMPCA$rotation
simpleCache("rsEnrichment_657", assignToVariable = "rsEnrichment")
simpleCache("rsEnrichmentTop10_657", assignToVariable = "rsEnrichmentTop10")
# one directory per PCA. This script generates plots for one PCA.
plotSubdir 
simpleCache(allMPCAString, assignToVariable = "allMPCA")
methylData = brcaMList$methylProp
# Sys.getenv("PLOTS") should be set

# CAPS variables should not be changed later in the script 
# (but R does not have immutable objects)
PCSTOANNOTATE = paste0("PC", 1:10)

#### plots that will be created and script specific parameters for them 
# "comparePCHeatmap"

# "methylAlongPC"
PCsToAnnotate_mAPC = PCSTOANNOTATE[1:5]

# "regionQuantileByPC"
PCsToAnnotate_rQBPC = c("PC1m4", "PC1p3", PCSTOANNOTATE)





##################################################################################
# comparePCHeatmap
# visualization of enrichment score results across PCs
# see if region set high in one PC is also high in others

### #script_specific
PCsToAnnotate = c("PC1m4", "PC1p3", PCSTOANNOTATE)


# number of plots = length(PCsToRankBy). one plot for each
# TODO: filter out low coverage region sets
# for rsEnrichment
comparePCHeatmap(rsEnrichment=rsEnrichment, 
                 PCsToRankBy=PCsToAnnotate, 
                 PCsToInclude=PCsToAnnotate,
                 fileName=paste0(Sys.getenv("PLOTS"), "rsEnrichHeatmap.pdf"))

# for rsEnrichmentTop10
comparePCHeatmap(rsEnrichment=rsEnrichmentTop10, 
                 PCsToRankBy=PCsToAnnotate, 
                 PCsToInclude=PCsToAnnotate,
                 fileName=paste0(Sys.getenv("PLOTS"), "rsEnrichHeatmapTop10Variable.pdf"))


##################################################################################
# methylAlongPC
# looking at methylation level data at individual cytosines ordered by PC 
# only looking at regions with high average loading scores
# still individual cytosine methylation

# one pdf for each PC given.
for (i in seq_along(PCsToAnnotate_mAPC)) {
    
    # top region sets for this PC
    rsInd = as.numeric(rsEnSortedInd[seq_along(topRSToPlotNum), PCsToAnnotate_mAPC[i], with=FALSE]) # original index
    
    grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir, "regionMethylHeatmaps", PCsToAnnotate_mAPC[i], ".pdf"), width = 11, height = 8.5 * topRSToPlotNum)
    
    # heatmap
    methylAlongPC(loadingMat=loadingMat, loadingThreshold=0.95, 
                  pcScores=allMPCA$x, 
                  coordinateDT=coordinateDT, 
                  methylData=methylData, 
                  GRList=GRList[rsInd], orderByPC=PCsToAnnotate_mAPC[i], 
                  topXRegions=50)
    dev.off()
}

###################################################################################
# regionQuantileByPC
# comparing loading scores/percentiles for individual regions among PCs
# need region sets and PCA loadings

# get top 10 region sets from each PC
topRSInd = unique(unlist(rsEnSortedInd[1:10, ]))

if (!dir.exists(paste0(Sys.getenv("PLOTS"), plotSubdir, "regionByPC/"))) {
    dir.create(paste0(Sys.getenv("PLOTS"), plotSubdir, "regionByPC/"))
}

grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir, 
                      "regionByPC/regionPercentileByPC", ".pdf"), 
               width = 11, height = 8.5 * length(topRSInd))

# ranking in terms of percentiles in case there were different distributions of loading scores for each PC
regionQuantileByPC(loadingMat=loadingMat, coordinateDT=coordinateDT, 
                   GRList=GRList[topRSInd], 
                   rsNames=paste0(rsEnrichment$rsName[topRSInd], " : ", rsEnrichment$rsDescription[topRSInd]), 
                   PCsToAnnotate=PCsToAnnotate_rQBPC)
dev.off()

##################################################################################
# seeing whether a subset of cytosines (ie a single region set) can
# recapitulate PC score from all cytosines
rsIndex = sort.int(rsEnrichment$PC2, index.return = TRUE, decreasing = TRUE)$ix[1:10]
regionSetList = GRList[rsIndex] 
regionSetName = paste0(rsEnrichment$rsName[rsIndex], " : ", rsEnrichment$rsDescription[rsIndex])
names(regionSetList) <-  regionSetName
subsetCorList = lapply(X = as.list(regionSetList), FUN = function(x) pcFromSubset(regionSet = x, 
                                                                                  pca = allMPCA, 
                                                                                  methylData = trainingMData, 
                                                                                  coordinateDT = brcaMList$coordinates, 
                                                                                  PCofInterest = paste0("PC", 1:10),
                                                                                  returnCor = TRUE))
i=5
Heatmap(subsetCorList[[i]], cluster_rows = FALSE, cluster_columns = FALSE, 
        column_title = names(subsetCorList)[i])
subsetCorList[[i]]
#cell_fun = function(j, i, x, y, width, height, fill) {
# grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
# })
