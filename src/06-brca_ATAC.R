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

source(ffCode("COCOA_paper/src/load_process_regions_brca.R"))
# reload my COCOA.R (from COCOA feat-atac branch)
source(ffCode("COCOA/R/COCOA.R"))
source(ffCode("COCOA/R/utility.R"))

##############################################################################

ryan_brca_count_matrix <- fread(ffData("tcga/ATACseq/TCGA-ATAC_BRCA_peaks_counts.tsv"), header=TRUE)

counts   <- as.matrix(ryan_brca_count_matrix[,2:75])
counts   <- counts[, order(colnames(counts))]
tcount   <- t(counts)
colnames(tcount) <- ryan_brca_count_matrix$sample

# Add metadata for each sample that has it
metadata <- fread(ffProc("COCOA_paper/analysis/atac/scores/brca/tcga_brca_metadata.csv"))
metadata <- metadata[!duplicated(metadata$subject_ID),]
metadata <- metadata[order(subject_ID),]

merged    <- as.data.frame(tcount)
merged$id <- rownames(tcount)
merged    <- merge(merged, metadata, by.x="id", by.y="subject_ID")

simpleCache(paste0("brcaATACPCA_", nrow(merged)), {
    pca       <- prcomp(as.matrix(merged[,2:215921]))
    pca
}, assignToVariable = "pca")

loadings  <- pca$rotation
scores    <- pca$x
rownames(scores) <- merged[,1]
brca_tcga_data   <- as.matrix(merged[,2:215921])

ggbiplot(pca, choices=c(1,2), groups=merged$ER_status, var.axes=FALSE, var.scale=1, obs.scale=1, ellipse=F, circle=F) + scale_color_discrete(name='') + theme(legend.direction='vertical', legend.position = 'right')
ggbiplot(pca, choices=c(1,6), groups=merged$ER_status, var.axes=FALSE, var.scale=1, obs.scale=1, ellipse=F, circle=F) + scale_color_discrete(name='') + theme(legend.direction='vertical', legend.position = 'right')


peaks           <- ryan_brca_count_matrix[, .(Chromosome, Start, End)]
colnames(peaks) <- c("chr","start","end")
pGR             <- makeGRangesFromDataFrame(peaks)

# Load methylation region sets
# Loads a GRList, rsName, and rsDescription object setls()
PCsToAnnotate = paste0("PC", 1:10)
system.time(
simpleCache(paste0("brcaATACrsScores_", nrow(merged)), {
    
    regionSetScoresTotal=runCOCOA(loadingMat=loadings, signalCoord=peaks, GRList=GRList, PCsToAnnotate=PCsToAnnotate, scoringMetric="regionMean", overlapMethod="total")
    regionSetScoresTotal$rsName = rsName
    regionSetScoresTotal$rsDescription = rsDescription
    regionSetScoresTotal
}, assignToVariable = "regionSetScoresTotal", recreate = TRUE)
)

system.time(
    simpleCache(paste0("brcaATACrsScoresWeightedMean_", nrow(merged)), {
        
        regionSetScoresTotal=runCOCOA(loadingMat=loadings, signalCoord=peaks, GRList=GRList, 
                                      PCsToAnnotate=PCsToAnnotate, 
                                      scoringMetric="regionMean", 
                                      overlapMethod="regionWeightedMean")
        regionSetScoresTotal$rsName = rsName
        regionSetScoresTotal$rsDescription = rsDescription
        regionSetScoresTotal
    }, assignToVariable = "regionSetScoresTotal", recreate = TRUE)
)



rsScores = regionSetScoresTotal
rssTotalComplete = regionSetScoresTotal

write_delim(head(rssTotalComplete[order(rssTotalComplete$PC1, decreasing=T),]), "TCGA-ATAC_BRCA_load-process-regions_TotalOM_headRegionSetScores.tsv",delim="\t")

rsScoreHeatmap(rssTotalComplete, PCsToAnnotate=paste0("PC", 1:4), rsNameCol = "rsName", orderByPC = "PC1", column_title = "Region sets ordered by score for PC1")

plotRSConcentration(rssTotalComplete, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="esr1|eralpha|eraa|gata3|foxa1|H3R17me")
plotRSConcentration(rssTotalComplete, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="ezh2|suz12")
plotRSConcentration(rssTotalComplete, colsToSearch = c("rsName", "rsDescription"), 
                    scoreColName = paste0("PC", 1:10), pattern="h3k9me3")


pdf("TCGA-ATAC_BRCA_regionSetScoresTotal_load-process-regions_rsConcentration_esr-eraa-eralpha.pdf", width=10, height=10)
plotRSConcentration(rssTotalComplete, colsToSearch = "rsName", pattern="esr|eraa|eralpha")
dev.off()


#########################################################################################################

# exploratory analysis
# setwd(ffProc("COCOA_paper/analysis/atac/scores/brca/tile_500-10_hg38"))
# countDT = fread("tcga_coverage.tab")

#######################################################################
# working locally

# aPCA = readRDS(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/TCGA-ATAC_BRCA_pca.Rds"))
aPCA = pca
# rsScores =  readRDS(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/TCGA-ATAC_BRCA_rssTotal.Rds"))
# load(ffProc("COCOA_paper/atac/load_this_john.RData"))
rsScores= regionSetScoresTotal
aMetadata = read.csv(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/tcga_brca_metadata.csv"))
load(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/brca_peak_pca_sample_names.RData")) # pcaNames
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
plot(as.matrix(pcScoreAnno[,c("PC1", "PC2")]))

#############################################################################

counts = read.table(file = "/home/jtl2hk/processed/COCOA_paper/TCGA-ATAC_BRCA_peaks_counts.tsv", sep = "\t", header = TRUE)
load(file="/home/jtl2hk/processed/COCOA_paper/RCache/topRSCombined.RData")
signalCoord = topRSCombined

counts[1:5, 1:5]
counts[1:5, 75:86]
countMetadata = counts[, 76:ncol(counts)]
coordinates = data.frame(chr = countMetadata$Chromosome, start = countMetadata$Start, end = countMetadata$End)
signalCoord = COCOA:::dtToGr(coordinates)
counts = counts[, 2:75]

atacPatients = colnames(counts)
atacPatients = sub(pattern = "....$", replacement = "", x = atacPatients)
atacPatients = gsub(pattern = ".", replacement = "-", x = atacPatients, fixed = TRUE)
setdiff(atacPatients, patientMetadata$subject_ID)
metadata = data.frame(subject_ID = atacPatients, stringsAsFactors = FALSE)
metadata = merge(metadata, patientMetadata, all.x=TRUE, sort=FALSE)

atacPCA = dimRedOnRS(regionSet = topRSCombined, methylData = counts, mCoord = signalCoord,  drMethod = "pca")
plot(atacPCA$x[, c(1,2)])
colnames(counts)
keepInd = atacPatients %in% patientMetadata$subject_ID
pcDF = atacPCA$x[keepInd, ]
pcDF = cbind(pcDF, metadata[1:70,])
colorClusterPlots(clusteredDF = pcDF, plotCols = c("PC1", "PC2"), colorByCols = c("ER_status", "PGR_status", "her2_status"))
plotPairwiseColPCs(pcaWithAnno = pcDF, 
                   pcsToPlot = paste0("PC", 1:10), 
                   colorByCols = c("ER_status"),#PGR_status", "her2_status"), 
                   plotDir = "/home/jtl2hk/processed/COCOA_paper/analysis/plots/atacPCA/")


atacTSNE = dimRedOnRS(regionSet = topRSCombined, methylData = counts, mCoord = signalCoord,  drMethod = "tsne", perplexity = 20)
plot(atacTSNE$Y)
tsneDF = atacTSNE$Y[keepInd, ]
colnames(tsneDF) = c("Dim1", "Dim2")
tsneDF = cbind(tsneDF, metadata[1:70, ])
colorClusterPlots(clusteredDF = tsneDF, plotCols = c("Dim1", "Dim2"), colorByCols = c("ER_status", "PGR_status", "her2_status"))

tsneDF[tsneDF$Dim1 > 1.25 & tsneDF$Dim2 < -3, ]


######### make "meta-region" loading profiles

# load region sets
# source(ffProjCode("load_process_regions_brca.R"))


topInd = rsRankingIndex(rsScores = rsScores, PCsToAnnotate = c("PC1", "PC2"))
topPC1Ind = topInd[, "PC1"][1:15]
topPC2Ind = topInd[, "PC2"][1:15]
uTopInd = unique(c(topPC1Ind, topPC2Ind))

# topRSList = GRList[uTopInd]
topRSList = lapply(X = GRList[uTopInd], FUN = function(x) resize(x = x, width = 14000, fix = "center"))
topRSNames = rsScores$rsName[uTopInd]


multiProfileP = makeMetaRegionPlots(loadingMat=pca$rotation, 
                                    signalCoord=peaks, GRList=topRSList, 
                                    rsNames=topRSNames, 
                                    PCsToAnnotate=PCsToAnnotate, binNum=21, overlapMethod="total") 
multiProfileP2 = makeMetaRegionPlots(loadingMat=pca$rotation, 
                                    signalCoord=peaks, GRList=topRSList, 
                                    rsNames=topRSNames, 
                                    PCsToAnnotate=PCsToAnnotate, binNum=21, overlapMethod="regionWeightedMean") 

inputID = paste0("brcaATAC", nrow(merged))
ggsave(filename = ffPlot(paste0(plotSubdir,
                         "/metaRegionLoadingProfiles", 
                         inputID, ".pdf")), plot = multiProfileP[[1]], device = "pdf", limitsize = FALSE)
ggsave(filename = ffPlot(paste0(plotSubdir, "/metaRegionLoadingProfilesWeightedMean", inputID, ".pdf")), plot = multiProfileP2[[1]], device = "pdf", limitsize = FALSE)

    