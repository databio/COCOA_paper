

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

#################################################################################
# place to save plots
plotSubdir = paste0(plotSubdir, "/")
if (!dir.exists(paste0(Sys.getenv("PLOTS"), plotSubdir))) {
    dir.create(paste0(Sys.getenv("PLOTS"), plotSubdir), recursive = TRUE)
}

source(paste0(Sys.getenv("CODE"), "aml_e3999/src/00-genericFunctions.R"))

inputID = addUnderscore(inputID, side = "left")

library(grid)
library(ggplot2)
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
                 fileName=paste0(Sys.getenv("PLOTS"), plotSubdir, 
                                 "rsEnrichHeatmap", inputID, ".pdf"))

##################################################################################
# methylAlongPC
# looking at methylation level data at individual cytosines ordered by PC 
# only looking at regions with high average loading scores
# still individual cytosine methylation

# one pdf for each PC given.
for (i in seq_along(PCsToAnnotate_mAPC)) {
    
    # top region sets for this PC
    rsInd = as.numeric(as.matrix(rsEnSortedInd[1:topRSToPlotNum, PCsToAnnotate_mAPC[i], with=FALSE])) # original index
    
    grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir, "regionMethylHeatmaps", PCsToAnnotate_mAPC[i], inputID, ".pdf"), width = 11, height = 8.5 * topRSToPlotNum)
    
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
                      "regionByPC/regionPercentileByPC", inputID, ".pdf"), 
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

.regionSetList = GRList[topRSInd_pcFSCH] 
regionSetName = paste0(rsEnrichment$rsName[topRSInd_pcFSCH], " : ", rsEnrichment$rsDescription[topRSInd_pcFSCH])
names(.regionSetList) <-  regionSetName
subsetCorList = lapply(X = as.list(.regionSetList), FUN = function(x) pcFromSubset(regionSet = x, 
                                                                                  pca = mPCA, 
                                                                                  methylData = methylData, 
                                                                                  coordinateDT = coordinateDT, 
                                                                                  PCofInterest = PCsToAnnotate_pcFSCHM,
                                                                                  returnCor = TRUE))
subsetCorMat = do.call(rbind, subsetCorList)
colnames(subsetCorMat) <- PCsToAnnotate_pcFSCHM

simpleCache("subsetCorMat", {subsetCorMat}, recreate=TRUE)

grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir, "subsetCorRSbyPC", inputID, ".pdf"), width = 8.5, 11)

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
.regionSetList = GRList[topRSInd_rsOLCP] 
pOL = percentCOverlap(coordGR = MIRA:::dtToGr(coordinateDT), 
                      GRList = .regionSetList) 

grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir, "topRSOverlap", inputID, ".pdf"), width = 25, height = 25)

Heatmap(matrix = pOL[[1]], cluster_rows = FALSE, cluster_columns = FALSE, 
        column_title = , cell_fun = function(j, i, x, y, width, height, fill, mat=pOL[[1]]) {
            grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()


####################################################################
# "meta region loading profiles" (mrLP)
# meta-region plots of region surrounding regions in region set
# check whether enrichment is specific to this region set by
# seeing if loading values have a spike in the center of these region sets
# compared to surrounding genome 
.regionSetList = GRList[topRSInd_mrLP]
.regionSetList = lapply(.regionSetList, resize, width = 14000, fix="center")

simpleCache(paste0("pcProf14k", inputID), {
    pcProf = pcEnrichmentProfile(loadingMat = mPCA$rotation, coordinateDT = coordinateDT,
                                 GRList=.regionSetList, PCsToAnnotate = PCsToAnnotate_mrLP,
                                 binNum = 21)
    # set names by reference
    #setattr(pcProf, "names", names(.regionSetList))
    pcProf
}, recreate = TRUE)
pcP = copy(get(paste0("pcProf14k", inputID)))

rsNames = paste0(rsEnrichment$rsName[topRSInd_mrLP], " : ", rsEnrichment$rsDescription[topRSInd_mrLP])



# average loading value from each PC to normalize so PCs can be compared with each other
avLoad = apply(X = loadingMat[, PCsToAnnotate_mrLP], MARGIN = 2, FUN = function(x) mean(abs(x)))

# normalize
# pcP = lapply(pcP, FUN = function(x) t(apply(X = x, MARGIN = 1, FUN = function(y) y - c(0, avLoad))))
pcP = lapply(pcP, FUN = function(x) x[, mapply(FUN = function(y, z) get(y) - z, y=PCsToAnnotate_mrLP, z = avLoad)])
pcP = lapply(pcP, FUN = function(x) data.table(regionGroupID=1:nrow(x), x))

# for the plot scale
maxVal = max(sapply(pcP, FUN = function(x) max(x[, get(PCsToAnnotate_mrLP)])))
minVal = min(sapply(pcP, FUN = function(x) min(x[, get(PCsToAnnotate_mrLP)])))

# convert to long format for plots
pcP = lapply(X = pcP, FUN = function(x) tidyr::gather(data = x, key = "PC", value="loading_value", PCsToAnnotate_mrLP))
pcP = lapply(X = pcP, as.data.table)
pcP = lapply(pcP, function(x) x[, PC := factor(PC, levels = PCsToAnnotate_mrLP)])

# xLabels = rep("", 21)
# xLabels[1] = "-14"
# xLabels[11] = "0"
# xLabels[21] = "14"

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
ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir,
                        "/metaRegionLoadingProfiles", 
                         inputID, ".pdf"), plot = multiProfileP, device = "pdf", limitsize = FALSE)

# check PPARG.bed, Jaspar motifs (had 18 rows instead of 21)
