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
# column was named sample but is actually region ID
colnames(tcount) <- ryan_brca_count_matrix$sample

# Add metadata for each sample that has it
metadata <- fread(ffProc("COCOA_paper/analysis/atac/scores/brca/tcga_brca_metadata.csv"))
metadata <- metadata[!duplicated(metadata$subject_ID),]
metadata <- metadata[order(subject_ID),]

merged    <- as.data.frame(tcount)
merged$id <- rownames(tcount)
merged    <- merge(merged, metadata, by.x="id", by.y="subject_ID")
pcaSampleNames <- merged$id

simpleCache(paste0("brcaATACPCA_", nrow(merged)), {
    pca       <- prcomp(as.matrix(merged[,2:215921]))
    row.names(pca$x) <- pcaSampleNames
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
