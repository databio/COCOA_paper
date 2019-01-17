# 



# project.init(codeRoot = paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/"), dataDir = paste0(Sys.getenv("PROCESSED"), "brca_PCA/"))
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
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

tiles = fread(paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/atac/tcga_coverage.tab"))
head(tiles)
nrow(tiles)
ncol(tiles)
tiles[, V1:= NULL]
tiles = apply(X = tiles, MARGIN = 2, FUN = as.numeric)
sampleTotals = colSums(tiles)
hist(colSums(tiles))

hist(tiles[tiles[,2] < 200,2], breaks = seq(0, 200, 2))
sum(tiles[tiles[,2] < 2,2])
max(tiles[,2])
hist(colMeans(tiles))


################################################################################3
# checking for outliers

rMean = apply(tiles, 1, mean) # mean for each region
hist(rMean[rMean < 200], breaks = seq(0, 200, 2))
medScore = median(rMean)
outlierNames = c(colnames(tiles)[apply(tiles[rMean < medScore,], 1, which.max)], 
                 colnames(tiles)[apply(tiles[rMean > medScore,], 1, which.min)])

#grDevices::pdf(file=paste0(Sys.getenv("PLOTS"), plotSubdir, "CpG_methyl_level_outliers_postQC2.pdf"))
plot(table(outlierNames), xlab="subject_ID", 
     ylab="Count of having most extreme ATAC level for a region", main="Samples with highest frequency of extreme methylation values")
dev.off()
