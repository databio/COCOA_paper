# below is based off code by Jason Smith

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(COCOA)
library("ComplexHeatmap")
library(ggbiplot)
library(readr)

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

loadGRList(genomeV = "hg38")
devtools::load_all(ffCode("COCOA/"))
loadBRCAatac(signalMat = TRUE, signalCoord = TRUE, 
             pcScores = TRUE, loadingMat = TRUE)

##############################################################################

rsScoreHeatmap(rssTotalComplete, PCsToAnnotate=paste0("PC", 1:4), rsNameCol = "rsName", orderByPC = "PC1", column_title = "Region sets ordered by score for PC1")

plotRSConcentration(rssTotalComplete, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="esr1|eralpha|eraa")
plotRSConcentration(rssTotalComplete, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="esr1|eralpha|eraa|gata3|foxa1|H3R17me")
plotRSConcentration(rssTotalComplete, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="ezh2|suz12")
plotRSConcentration(rssTotalComplete, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="h3k9me3")


pdf("TCGA-ATAC_BRCA_regionSetScoresTotal_load-process-regions_rsConcentration_esr-eraa-eralpha.pdf", width=10, height=10)
plotRSConcentration(rssTotalComplete, scoreColName = "PC1", colsToSearch = "rsName", pattern="esr|eraa|eralpha")
dev.off()

erPC1Hist = plotRSConcentration(rssTotalComplete, 
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


topInd = rsRankingIndex(rsScores = rsScores, signalCol = c("PC1", "PC2"))
topPC1Ind = topInd[, "PC1"][1:15]
topPC2Ind = topInd[, "PC2"][1:15]
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

    