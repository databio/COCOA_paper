
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
plotSubdir = "23-brcaATAC/"

if(!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}


# cache data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))

regionCovCutoff = 100

##########

loadBRCAatac(signalMat = TRUE, signalCoord = TRUE, 
             pcScores = TRUE, loadingMat = TRUE)

# parameters that should have been loaded by previous script
# (load them now if I'm just running this script by itself)
if (!exists("nPerm")) {
    nPerm = 300
}
if (!exists("dataID")) {
    dataID =paste0("brcaATAC", ncol(signalMat))
    # dataID = "brcaATACPCAPerm"
}
if (!exists("variationMetric")) {
    variationMetric = "cor"
}

.analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)
inputID = paste0("_", nPerm, "Perm_", variationMetric, 
                                     "_", "brcaATAC", ncol(signalMat))

simpleCache(paste0("rsScores", paste0("_", "brcaATAC", ncol(signalMat)), 
                   "_", variationMetric), assignToVariable = "rsScores", reload = TRUE)

# simpleCache(paste0("pRankedScores", .analysisID), assignToVariable = "rsScores", reload = TRUE)
# rsEnSortedInd = rsRankingIndex(rsScores = rsScores, 
#                                signalCol = list(paste0(paste0("PC", 1:10), "_PValGroup"), 
#                                                 paste0("PC", 1:10)), 
#                                decreasing = c(TRUE, TRUE), newColName = paste0("PC", 1:10))
simpleCache(paste0("rsPermScores", .analysisID), assignToVariable = "rsPermScores")
signalCol = paste0("PC", 1:10)

loadGRList(genomeV = "hg38")

# screen out region sets that have less than 100 regions that overlap the data
if (length(GRList) != nrow(rsScores)) {
    stop("GRList and rsScores do not match")
}
lowCovRS <- rsScores$regionSetCoverage < 100
rsScores <- rsScores[!lowCovRS, ]
GRList <- GRList[!lowCovRS]
rsName <- rsName[!lowCovRS]
rsDescription <- rsDescription[!lowCovRS]

rsEnSortedInd = rsRankingIndex(rsScores = rsScores, 
                               signalCol = paste0("PC", 1:10), 
                               decreasing = TRUE)

##############################################################################
# panel A
########## make PCA plot
aMetadata = read.csv(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/tcga_brca_metadata.csv"))
pcaNames = row.names(aPCA$x)
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
# panel B
# load(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/brca_peak_pca_sample_names.RData")) # pcaNames

pcAnnoScoreDist = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = "PC1", 
                  pattern = c("esr|eralpha", "foxa1|gata3|H3R17me2"), 
                  patternName = c("ER", "ER-related")) + 
                 theme(text = element_text(size=15), legend.position = c(0, 0.2)) 
                  # coord_fixed(ratio = 10)
pcAnnoScoreDist 
ggsave(filename = ffPlot(paste0(plotSubdir, "pc1AnnoScoreDistERRelated.svg")), 
       plot = pcAnnoScoreDist, device = "svg", height = 100, width = 100, units = "mm")
# inkscape uses 90 dpi instead of 72
# inkscape total plot size dimensions are 0.8 times ggplot dimensions
pcAnnoScoreDist2 = plotAnnoScoreDist2(rsScores = rsScores, colsToPlot = "PC1", pattern = c("esr|eralpha", "foxa1|gata3|H3R17me2"), 
                   patternName = c("ER", "ER-related"))
pcAnnoScoreDist2 
                                      #pattern = "esr|eralpha|foxa1|gata3|H3R17me2", patternName = "ER-related")
ggsave(filename = ffPlot(paste0(plotSubdir, "pc1AnnoScoreDist2ERRelated.svg")), 
       plot = pcAnnoScoreDist2, device = "svg",  height = 150, width = 150, units = "mm")


################# 
# panel C 
# PC2, immune-related
# reviews of hematopoietic transcription factors (not exhaustive obviously):
# one general review, one myeloid review, and one lymphoid review
# myeloid: https://www.nature.com/articles/nri2024
# (table 1) RUNX1, SCL/TAL1, PU.1, CEBPa, IRF8, GFI1, CEBPe 

# abbrev cebpa
hemaTFs = c("RUNX1", "SCL|TAL1", "PU.1|PU1|SPI1", 
            "CEBPA", "IRF8", "GFI1", "CEBPE")

# general: http://www.jbc.org/content/270/10/4955.short
# (table)
# TCF3 (E2A), KLF1 (EKLF), GATA1, GATA2, Ikaros, c-MYB, p45 NF-E2, PAX5, PU.1, RBTN2, SCL/TAL1
hemaTFs = c(hemaTFs, c("TCF3", "KLF1", "GATA1", "GATA2", "Ikaros|IKZF1", "CMYB", "NFE2"))

# lymphoid: https://doi.org/10.1182/blood-2014-12-575688
# (figure 1. table 1) IKZF1, TCF3, EBF1, PAX5, FOXO1, ID2, GATA3
hemaTFs = c(hemaTFs, c("TCF3", "EBF1", "PAX5", "FOXO1", "ID2", "GATA3"))

hemaTFs = unique(hemaTFs)

hemaPattern = paste0(hemaTFs, collapse = "|")
plotAnnoScoreDist(rsScores = rsScores, colsToPlot = "PC2", pattern = hemaPattern, patternName = "Hematopoietic TFs")
pc2AnnoScoreDist = plotAnnoScoreDist2(rsScores = rsScores, colsToPlot = "PC2", 
                   pattern = hemaPattern, patternName = "Hematopoietic TFs")
ggsave(ffPlot(paste0(plotSubdir, "pc2HemaATAC.svg")), 
       plot = pc2AnnoScoreDist, device = "svg")


######### make "meta-region" loading profiles
# panel D
atacCor = cor(x = t(signalMat), y = pcScores[, signalCol])
all(colnames(signalMat) == row.names(pcScores))

# load region sets
# source(ffProjCode("load_process_regions_brca.R"))

topPC1Ind = rsEnSortedInd[, "PC1"][1:15]
topPC2Ind = rsEnSortedInd[, "PC2"][1:15]
uTopInd = unique(c(topPC1Ind, topPC2Ind))

# topRSList = GRList[uTopInd]
topRSList = lapply(X = GRList[uTopInd], FUN = function(x) resize(x = x, width = 14000, fix = "center"))
topRSNames = rsScores$rsName[uTopInd]


multiProfileP = makeMetaRegionPlots(signal=atacCor, 
                                    signalCoord=signalCoord, GRList=topRSList, 
                                    rsNames=topRSNames, 
                                    signalCol=signalCol, binNum=21, aggrMethod ="simpleMean") 
multiProfileP2 = makeMetaRegionPlots(signal=atacCor, 
                                    signalCoord=signalCoord, GRList=topRSList, 
                                    rsNames=topRSNames, 
                                    signalCol=signalCol, binNum=21, aggrMethod ="proportionWeightedMean") 

ggsave(filename = ffPlot(paste0(plotSubdir,
                         "/metaRegionLoadingProfiles", 
                         inputID, ".pdf")), plot = multiProfileP[[1]], device = "pdf", limitsize = FALSE)
ggsave(filename = ffPlot(paste0(plotSubdir, "/metaRegionLoadingProfilesWeightedMean", inputID, ".pdf")), 
       plot = multiProfileP2[[1]], device = "pdf", limitsize = FALSE)

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

