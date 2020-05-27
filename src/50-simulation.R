## generate simulated DNA methylation data to compare COCOA and LOLA

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(curatedTCGAData)
library(TCGAutils)

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))

genomeV = "hg19"

set.seed(1234)


##################################
# start from single healthy sample
loadGRList(genomeV = genomeV)

loadTCGAMethylation(cancerID = "KIRC")
methylMat = methylList$methylProp
signalCoord = methylList$coordinates

sampleType = substr(colnames(methylMat), start = 14, stop = 15)
# 01 is primary solid tumor, 11 is solid normal tissue, 05 is new primary tumor
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
normalSampleInd = (sampleType == "11")
tumorSampleInd = (sampleType == "01") # exclude the extra sample 05
methylMat = methylMat[, normalSampleInd]
# now I only need patient ID
colnames(methylMat) = substr(colnames(methylMat), start = 1, stop = 12)

# average of all normal samples (160 samples)
healthy = apply(X = methylMat, 1, mean)

# fake tumor profile, only changing DNA methylation in one region set
tumor = healthy

# find CpGs that overlap region set
signalCoord = COCOA:::dtToGr(signalCoord)
myGR = GRList[[c("Human_MDA-MB-231-Cells_ESR1,-DBDmut_E2-45min_Katzenellenbogen.bed")]] # 356 
fo = findOverlaps(query = myGR, subject = signalCoord)
olCGInd = subjectHits(fo)
tumor[olCGInd] = 1
healthy[olCGInd] = 0



# H1hesc_WE.bed




# set healthy sample to opposite of tumor in this region set
    

###################################
sampleN = 20
nCG = 10000
trueProp = 0.05
cgProps = rep(-1, nCG)
for (i in 1:nCG) {
  cgProps[i] = mean(sample(c(0, 1), size = sampleN, replace = TRUE, prob = c(1-trueProp, trueProp)))
}
hist(cgProps, breaks = seq(0, 1, by=0.05))
