


# project.init(codeRoot = paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/"), dataDir = paste0(Sys.getenv("PROCESSED"), "brca_PCA/"))
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(ggplot2)
# library(fastICA)

# 
setwd(paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/plots/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)


# DNA methylation data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))

#############################################################################

counts = read.table(file = "/home/jtl2hk/processed/brca_PCA/TCGA-ATAC_BRCA_peaks_counts.tsv", sep = "\t", header = TRUE)
load(file="/home/jtl2hk/processed/brca_PCA/RCache/topRSCombined.RData")
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
                   plotDir = "/home/jtl2hk/processed/brca_PCA/analysis/plots/atacPCA/")


atacTSNE = dimRedOnRS(regionSet = topRSCombined, methylData = counts, mCoord = signalCoord,  drMethod = "tsne", perplexity = 20)
plot(atacTSNE$Y)
tsneDF = atacTSNE$Y[keepInd, ]
colnames(tsneDF) = c("Dim1", "Dim2")
tsneDF = cbind(tsneDF, metadata[1:70, ])
colorClusterPlots(clusteredDF = tsneDF, plotCols = c("Dim1", "Dim2"), colorByCols = c("ER_status", "PGR_status", "her2_status"))

tsneDF[tsneDF$Dim1 > 1.25 & tsneDF$Dim2 < -3, ]

################################################################################
#                                                                              #
#                                   START                                      #
#                                                                              #
################################################################################
# below is based off code by Jason Smith


library(LOLA)
library(COCOA)
library(data.table)
library(ggplot2)
library("ComplexHeatmap")
library(ggbiplot)

cocoa_dir <- "/scratch/jps3dp/tools/databio/COCOA/R/" # feat-atac branch
data_dir  <- "/sfs/lustre/allocations/shefflab/processed/COCOA_paper/analysis/atac/"
tcga_dir  <- "/scores/brca/tcga_brca_peaks-log2counts-dedup/"

source(paste0(cocoa_dir, "COCOA.R"))
source(paste0(cocoa_dir, "utility.R"))

ryan_brca_count_matrix <- fread("/sfs/lustre/allocations/shefflab/data/tcga/ATACseq/TCGA-ATAC_BRCA_peaks_counts.tsv", header=TRUE)

counts   <- as.matrix(ryan_brca_count_matrix[,2:75])
counts   <- counts[, order(colnames(counts))]
tcount   <- t(counts)
colnames(tcount) <- ryan_brca_count_matrix$sample

# Add metadata for each sample that has it
metadata <- fread("/sfs/lustre/allocations/shefflab/processed/COCOA_paper/analysis/atac/scores/brca/tcga_brca_metadata.csv")
metadata <- metadata[!duplicated(metadata$subject_ID),]
metadata <- metadata[order(subject_ID),]

merged    <- as.data.frame(tcount)
merged$id <- rownames(tcount)
merged    <- merge(merged, metadata, by.x="id", by.y="subject_ID")

pca       <- prcomp(as.matrix(merged[,2:215921]))
loadings  <- pca$rotation
scores    <- pca$x
rownames(scores) <- merged[,1]
brca_tcga_data   <- as.matrix(merged[,2:215921])

ggbiplot(pca, choices=c(1,2), groups=merged$ER_status, var.axes=FALSE, var.scale=1, obs.scale=1, ellipse=F, circle=F) + scale_color_discrete(name='') + theme(legend.direction='vertical', legend.position = 'right')

peaks           <- ryan_brca_count_matrix[,76:78]
colnames(peaks) <- c("chr","start","end")
pGR             <- makeGRangesFromDataFrame(peaks)

# Load methylation region sets
# Loads a GRList, rsName, and rsDescription object setls()

source("/sfs/lustre/scratch/jps3dp/tools/databio/COCOA_paper/src/load_process_regions_brca.R")
# reload my COCOA.R (from COCOA feat-atac branch)
source("/sfs/lustre/scratch/jps3dp/tools/databio/COCOA/R/COCOA.R")
source("/sfs/lustre/scratch/jps3dp/tools/databio/COCOA/R/utility.R")

system.time(regionSetScoresTotal <- runCOCOA(loadingMat=loadings, signalCoord=peaks, GRList=GRList, PCsToAnnotate=PCsToAnnotate, scoringMetric="regionMean", overlapMethod="total"))

library(readr)
write_delim(head(rssTotalComplete[order(rssTotalComplete$PC1, decreasing=T),]), "TCGA-ATAC_BRCA_load-process-regions_TotalOM_headRegionSetScores.tsv",delim="\t")

rsScoreHeatmap(rssTotalComplete, PCsToAnnotate=paste0("PC", 1:4), rsNameCol = "regionSetName", orderByPC = "PC1", column_title = "Region sets ordered by score for PC1")

plotRSConcentration <- function(rsScores, scoreColName="PC1", 
                                colsToSearch = c("rsName", "rsDescription"), 
                                pattern, percent = FALSE) {
    # breaks
    
    rsRankInd = rsRankingIndex(rsScores=rsScores, PCsToAnnotate=scoreColName)
    
    
    rsInd = rep(FALSE, nrow(rsScores))
    for (i in seq_along(colsToSearch)) {
        rsInd = rsInd | grepl(pattern = pattern, x = rsScores[, colsToSearch[i]], ignore.case = TRUE)
    }
    
    rsScores$ofInterest = rsInd
    ofInterestDF = as.data.frame(rsInd[as.matrix(rsRankInd)])
    colnames(ofInterestDF) <- colnames(rsRankInd)
    ofInterestDF$rsRank = 1:nrow(ofInterestDF)
    categoryDistPlot = ggplot(ofInterestDF, aes(x=rsRank, weight=get(scoreColName))) + 
        geom_histogram() + theme_classic()#+ facet_wrap(~get(scoreColName))
    return(categoryDistPlot)
    
}

plotRSConcentration(rssTotalComplete, colsToSearch = "regionSetName", pattern="esr")

pdf("TCGA-ATAC_BRCA_regionSetScoresTotal_load-process-regions_rsConcentration_esr-eraa-eralpha.pdf", width=10, height=10)
plotRSConcentration(rssTotalComplete, colsToSearch = "regionSetName", pattern="esr|eraa|eralpha")
dev.off()

################################################################################
#                                                                              #
#                                    END                                       #
#                                                                              #
################################################################################