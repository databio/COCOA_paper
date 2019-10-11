
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(COCOA)
library("ComplexHeatmap")
library(ggbiplot)
library(readr)
devtools::load_all(ffCode("COCOA/"))

# 
setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)
plotSubdir = "06-brcaATAC/"

if(!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}


# DNA methylation data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))

##########

cocoa_dir <- ffCode("COCOA/R/") # feat-atac branch
data_dir  <- ffProc("COCOA_paper/analysis/atac/")
tcga_dir  <- "/scores/brca/tcga_brca_peaks-log2counts-dedup/"



loadBRCAatac(signalMat = TRUE, signalCoord = TRUE, 
             pcScores = TRUE, loadingMat = TRUE)
# nPerm = 250
# dataID =paste0("brcaATAC", ncol(signalMat))
# dataID = "brcaATACPCAPerm"
# variationMetric = "cor"
# .analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)
inputID = paste0("_", nPerm, "Perm_", variationMetric, 
                                     "_", "brcaATAC", ncol(signalMat))
simpleCache(paste0("pRankedScores", .analysisID), assignToVariable = "rsScores", reload = TRUE)

rsEnSortedInd = rsRankingIndex(rsScores = rsScores, 
                               signalCol = list(paste0(paste0("PC", 1:10), "_PValGroup"), 
                                                paste0("PC", 1:10)), 
                               decreasing = c(TRUE, TRUE), newColName = paste0("PC", 1:10))
#dataID =paste0("brcaATAC", ncol(signalMat))
.analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)
simpleCache(paste0("rsPermScores", .analysisID), assignToVariable = "rsPermScores")


loadGRList(genomeV = "hg38")
##############################################################################

rsScoreHeatmap(rsScores, signalCol=paste0("PC", 1:4), rsNameCol = "rsName", orderByCol = "PC1", column_title = "Region sets ordered by score for PC1")

plotRSConcentration(rsScores, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="esr1|eralpha|eraa")
plotRSConcentration(rsScores, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="esr1|eralpha|eraa|gata3|foxa1|H3R17me")
plotRSConcentration(rsScores, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="ezh2|suz12")
plotRSConcentration(rsScores, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="h3k9me3")


pdf("TCGA-ATAC_BRCA_regionSetScoresTotal_load-process-regions_rsConcentration_esr-eraa-eralpha.pdf", width=10, height=10)
plotRSConcentration(rsScores, scoreColName = "PC1", colsToSearch = "rsName", pattern="esr|eraa|eralpha")
dev.off()

erPC1Hist = plotRSConcentration(rsScores, 
                    scoreColName = "PC1", 
                    colsToSearch = c("rsName"), 
                    pattern="esr|eraa|eralpha") + theme_get() + theme(axis.text.y= element_text(size=15),
                                                        axis.text.x = element_text(size=15),
                                                        axis.title.x = element_text(size=20),
                                                        axis.title.y = element_text(size=20), 
                                                        plot.title = element_text(size=20, hjust = 0.5))

ggsave(filename = ffPlot(paste0(plotSubdir, "rsConcentrationERPC1_brcaATAC_", nrow(merged), ".svg")), 
       plot = erPC1Hist, device = "svg") 

#########################################################################################################
# visualize region set distribution

aPCA = pca
rsScores= regionSetScoresTotal
aMetadata = read.csv(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/tcga_brca_metadata.csv"))
# load(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/brca_peak_pca_sample_names.RData")) # pcaNames
dim(aMetadata)
length(unique(aMetadata$subject_ID))
table(aMetadata$ER_status)
View(rsScores[order(rsScores$PC2, decreasing = TRUE), ])

plotRSConcentration(rsScores = rsScores, paste0("PC", 1:4), colsToSearch = "rsName", pattern = "esr|eralpha|foxa1|gata3|H3R17me2")
plotRSConcentration(rsScores = rsScores, "PC1", colsToSearch = "rsName", pattern = "ezh2|suz12")
plotRSConcentration(rsScores = rsScores, scoreColName = "PC2", colsToSearch = "rsName", pattern = "ezh2|suz12|h3k27me")
plotRSConcentration(rsScores = rsScores, scoreColName = "PC3", colsToSearch = "rsName", pattern = "ezh2|suz12")
plotRSConcentration(rsScores = rsScores, scoreColName = paste0("PC", 1:4), colsToSearch = "rsName", pattern = "ezh2|suz12")
plotRSConcentration(rsScores = rsScores, scoreColName = paste0("PC", 1:4), colsToSearch = "rsName", pattern = "runx|tal1|gata|lmo2|PU1|lyl1|FLI1|evi1")
pc1ERATAC = plotRSConcentration(rsScores = rsScores, "PC1", colsToSearch = "rsName", pattern = "esr|eralpha|foxa1|gata3|H3R17me2")
pc1ERATAC + theme(axis.text = element_text(colour = "black", size = 15), axis.ticks = element_line(colour = "black"))
pc1ERATAC
ggsave(ffPlot(paste0(plotSubdir, "pc1ERATAC.svg")), plot = pc1ERATAC, device = "svg")

################# PC2, immune-related
# reviews of hematopoietic transcription factors (not exhaustive obviously):
# myeloid: https://www.nature.com/articles/nri2024
# (table 1) RUNX1, SCL/TAL, PU.1, CEBPa, IRF8, GFI1, CEBPe 
# need to finish

# lymphoid: https://doi.org/10.1182/blood-2014-12-575688
# 

# general: http://www.jbc.org/content/270/10/4955.short
# (table)
# TCF3 (E2A), KLF1 (EKLF), GATA1, GATA2, Ikaros, c-MYB, p45 NF-E2, PAX5, PU.1, RBTN2, SCL/TAL1



########## make PCA plot

pcScore = data.frame(aPCA$x, pcaNames)
pcScoreAnno= merge(pcScore, aMetadata, by.x = "pcaNames", by.y= "subject_ID", all.x=TRUE)
colorClusterPlots(pcScoreAnno, plotCols = paste0("PC", c(1,2)), colorByCols = "ER_status")
ggplot(data = pcScoreAnno, mapping = aes(x = PC1, y= PC2)) + geom_point(aes(col=ER_status), size = 4, alpha=0.5) + theme_classic()
pcScoreAnno$ER_status[pcScoreAnno$ER_status == ""] = NA
aPCAPlot = ggplot(data = pcScoreAnno, mapping = aes(x = PC1, y= PC2)) + geom_point(aes(col=ER_status), size = 4, alpha=0.5) +
    theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20)) + coord_fixed()
aPCAPlot
ggsave(filename = ffPlot(paste0(plotSubdir, "pc1_2_BRCA_ATAC.svg")), plot = aPCAPlot, device = "svg")

plot(as.matrix(pcScoreAnno[,c("PC1", "PC2")]))

#############################################################################
######### make "meta-region" loading profiles

# load region sets
# source(ffProjCode("load_process_regions_brca.R"))


topPC1Ind = rsEnSortedInd[, "PC1"][1:15]
topPC2Ind = rsEnSortedInd[, "PC2"][1:15]
uTopInd = unique(c(topPC1Ind, topPC2Ind))

# topRSList = GRList[uTopInd]
topRSList = lapply(X = GRList[uTopInd], FUN = function(x) resize(x = x, width = 14000, fix = "center"))
topRSNames = rsScores$rsName[uTopInd]


multiProfileP = makeMetaRegionPlots(signal=pca$rotation, 
                                    signalCoord=peaks, GRList=topRSList, 
                                    rsNames=topRSNames, 
                                    signalCol=PCsToAnnotate, binNum=21, aggrMethod ="simpleMean") 
multiProfileP2 = makeMetaRegionPlots(signal=pca$rotation, 
                                    signalCoord=peaks, GRList=topRSList, 
                                    rsNames=topRSNames, 
                                    signalCol=PCsToAnnotate, binNum=21, aggrMethod ="proportionWeightedMean") 

inputID = paste0("brcaATAC", nrow(merged))
ggsave(filename = ffPlot(paste0(plotSubdir,
                         "/metaRegionLoadingProfiles", 
                         inputID, ".pdf")), plot = multiProfileP[[1]], device = "pdf", limitsize = FALSE)
ggsave(filename = ffPlot(paste0(plotSubdir, "/metaRegionLoadingProfilesWeightedMean", inputID, ".pdf")), plot = multiProfileP2[[1]], device = "pdf", limitsize = FALSE)

############################################################################

realRSScores = rsScores
# include null distributions
permResultsMat = do.call(cbind, lapply(rsPermScores, function(x) x$PC1))

# make a ROC curve plot for EZH2/Suz12
# pred = realRSScores$cancerStage / max(realRSScores$cancerStage, na.rm = TRUE)
pred = cbind(realRSScores$PC1, permResultsMat)
rocPreds = ROCR::prediction(pred, labels = matrix(grepl(pattern = "esr1|eralpha", 
                                                        x = realRSScores$rsName, 
                                                        ignore.case = TRUE), nrow=nrow(pred), ncol=ncol(pred)))
testAUC = ROCR::performance(rocPreds, measure="auc")@y.values[[1]]
allTestAUC = unlist(ROCR::performance(rocPreds, measure="auc")@y.values)
testAUC
perf = ROCR::performance(rocPreds, measure = "tpr", x.measure = "fpr")
plot(perf)
title(main= "ROC curve for estrogen receptor")

# add real score
pred = realRSScores$PC1
rocPreds = ROCR::prediction(pred, labels =grepl(pattern = "esr1|eralpha", 
                                                x = realRSScores$rsName, 
                                                ignore.case = TRUE) )
perf = ROCR::performance(rocPreds, measure = "tpr", x.measure = "fpr")
plot(perf, add=TRUE, col="red", lwd=3)

sort(allTestAUC)
hist(allTestAUC)
dev.off()

plot(perf, col="red", lwd=3)
# shuffling ER label of region sets
erLabel = grepl(pattern = "esr1|eralpha", 
                x = realRSScores$rsName, 
                ignore.case = TRUE)
for (i in 1:1000) {
    rocPreds = ROCR::prediction(pred, labels =sample(x = erLabel, size = length(erLabel), replace=FALSE))
    perf = ROCR::performance(rocPreds, measure = "tpr", x.measure = "fpr")
    plot(perf, add=TRUE)
}

################################################################
# add real score
pred = rsPermScores[[which.max(allTestAUC) - 1]]$PC1
rocPreds = ROCR::prediction(pred, labels =grepl(pattern = "esr1|eralpha", 
                                                x = realRSScores$rsName, 
                                                ignore.case = TRUE) )
perf = ROCR::performance(rocPreds, measure = "tpr", x.measure = "fpr")
plot(perf, col="red", lwd=3)

sort(allTestAUC)
hist(allTestAUC)
dev.off()

plot(perf, col="red", lwd=3)
# shuffling ER label of region sets
erLabel = grepl(pattern = "esr1|eralpha", 
                x = realRSScores$rsName, 
                ignore.case = TRUE)
for (i in 1:1000) {
    rocPreds = ROCR::prediction(pred, labels =sample(x = erLabel, size = length(erLabel), replace=FALSE))
    perf = ROCR::performance(rocPreds, measure = "tpr", x.measure = "fpr")
    plot(perf, add=TRUE)
}

# mPCA

############################################################################

# # testing if macrophage marker is correlated with PC2
# sharedSamples = intersect(row.names(pcScores), colnames(exprMat))
# cor.test(x = pcScores[sharedSamples, "PC2"], y = exprMat["CD68", sharedSamples]) # general marker
# cor.test(x = pcScores[sharedSamples, "PC2"], y = exprMat["CD163L1", sharedSamples], method="spearman") # 
# grep("CD1", row.names(exprMat), value=TRUE)

############################################################################
rsEnrichment = rsScores
coordinateDT = COCOA::grToDt(signalCoord)
mPCA = list()
mPCA$rotation = loadingMat
mPCA$x = pcScores
mPCA$center = rowMeans(methylData)
methylData = signalMat

PCSTOANNOTATE = paste0("PC", 1:10)


### plots that will be created and script specific parameters for them 
# "comparePCHeatmap"
PCsToAnnotate_cPCH = PCSTOANNOTATE
# "methylAlongPC"
topRSToPlotNum = 15
PCsToAnnotate_mAPC = PCSTOANNOTATE[1:10]
# "regionQuantileByPC"
PCsToAnnotate_rQBPC = PCSTOANNOTATE
topRSInd_rQBPC = unique(unlist(rsEnSortedInd[1:15, ])) # get top region sets from each PC
# pcFromSubset Correlation Heatmap
PCsToAnnotate_pcFSCH = PCSTOANNOTATE
topRSInd_pcFSCH = unique(unlist(rsEnSortedInd[1:15, ])) # get top region sets from each 
## "region set Overlapping Cytosine Proportion" (rsOLCP)
## proportion of cytosines from region set that are shared with other region set
topRSInd_rsOLCP = unique(unlist(rsEnSortedInd[1:10, ]))
## "meta region loading profiles" (mrLP)
topRSInd_mrLP = unique(unlist(rsEnSortedInd[1:10, ]))
PCsToAnnotate_mrLP = PCSTOANNOTATE

# the pipeline
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/COCOA_vis_pipeline.R")) 

