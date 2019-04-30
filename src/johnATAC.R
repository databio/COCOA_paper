# project.init(codeRoot = paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/"), dataDir = paste0(Sys.getenv("PROCESSED"), "brca_PCA/"))
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
#patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)
plotSubdir = "jATAC/"

#########################################################################################################

library(data.table)
setwd("/sfs/lustre/allocations/shefflab/processed/COCOA_paper/analysis/atac/scores/brca/tile_500-10_hg38")

countDT = fread("tcga_coverage.tab")

dim(countDT)
head(countDT)
# row.names(countDT) = countDT[, 1]
countDT[, V1 := NULL]

# region sums

# sample sums
cSum = colSums(countDT)
rSum = rowSums(countDT)
sum(rSum == 0)

max(countDT)

#######################################################################
# working locally

aPCA = readRDS(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/TCGA-ATAC_BRCA_pca.Rds"))
rsScores =  readRDS(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/TCGA-ATAC_BRCA_rssTotal.Rds"))
aMetadata = read.csv(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/tcga_brca_metadata.csv"))
dim(aMetadata)
length(unique(aMetadata$subject_ID))
table(aMetadata$ER_status)

plotRSConcentration(rsScores = rsScores, "PC1", colsToSearch = "regionSetName", pattern = "esr|eralpha|foxa1|gata3|H3R17me2")

plotRSConcentration(rsScores = rsScores, scoreColName = "PC2", colsToSearch = "regionSetName", pattern = "esr|eralpha|foxa1|gata3|H3R17me2")

plotRSConcentration(rsScores = rsScores, scoreColName = "PC4", colsToSearch = "regionSetName", pattern = "esr|eralpha|foxa1|gata3|H3R17me2")

########## make PCA plot
load("/home/jtl2hk/processed/COCOA_paper/atac/brca_peak_pca_sample_names.RData")
