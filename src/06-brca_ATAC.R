


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

