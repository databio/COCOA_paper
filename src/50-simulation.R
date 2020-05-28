## generate simulated DNA methylation data to compare COCOA and LOLA

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(curatedTCGAData)
library(TCGAutils)
library(bumphunter)
library(LOLA)

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))

genomeV = "hg19"
dataID = "KIRC_simulated"

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
# # H1hesc_WE.bed
fo = findOverlaps(query = myGR, subject = signalCoord)
olCGInd = subjectHits(fo)
tumor[olCGInd] = 1
# set healthy sample to opposite of tumor in this region set
healthy[olCGInd] = 0

# make samples with varying proportion healthy vs "tumor"
# name says what percent tumor
percTumor = seq(10, 100, 10)
percHealth = 100 - percTumor
mixed = ((tumor %*% t(percTumor)) + (healthy %*% t(percHealth))) / 100
colnames(mixed) = paste0("tumor", percTumor)

both = matrix(data = healthy, nrow = length(healthy), ncol=10)
colnames(both) = paste0("healthy", 1:10)
# 10 healthy, 10 tumor
both = cbind(both, mixed) 
bothPCA = prcomp(t(both))$x
plot(bothPCA[, c(1, 2)])

smallChange = runCOCOA(genomicSignal = both, signalCoord = signalCoord, GRList = GRList, 
         signalCol = c("PC1", "PC2"), targetVar = bothPCA, 
         variationMetric = "cov", scoringMetric = "regionMean", 
         absVal = TRUE, centerGenomicSignal = TRUE, centerTargetVar = TRUE)
smallChange = cbind(smallChange, rsName)

View(arrange(smallChange, desc(PC1)))

###############
# get DMRs with bumphunter
sampleStatus = as.numeric(grepl(pattern = "tumor", 
                               x = colnames(both), ignore.case = TRUE))

# Function gets DMRs with bumphunter then runs LOLA on those
# @param signalCoord GRanges.
# @param GRList GrangesList for LOLA
# @param rsAnno data.frame. One row for each region set in GRList.
# @param targetVar numeric. Can be discrete or continuous. Used with bumphunter. 
dmrLOLA = function(genomicSignal, signalCoord, GRList, rsAnno, targetVar, dataID="") {
    
    both = genomicSignal
    sampleStatus = targetVar
    
    designMat = model.matrix(object = ~ sampleStatus)
    tumorDMR =bumphunter(object = both, design=designMat, chr=COCOA:::grToDt(signalCoord)$chr, 
                         pos=COCOA:::grToDt(signalCoord)$start, coef=2, cutoff = 0.1, B=1000)# , type="Beta")
    
    simpleCache(paste0("tumorDMR_", dataID), {
        tumorDMR
    }, assignToVariable = "tumorDMR")
    View(tumorDMR$table)
    
    # LOLA
    dmrGR = COCOA:::dtToGr(tumorDMR$table[tumorDMR$table$fwer <=0.05, c("chr", "start", "end")])
    uniGR = resize(x = signalCoord, width = 2000, fix = "center")
    uniGR = reduce(uniGR)
    names(GRList) = NULL
    regionSetDB = list()
    regionSetDB$regionAnno = rsAnno
    regionSetDB$regionGRL = GRList
    lResults = LOLA::runLOLA(userSets = dmrGR, userUniverse = uniGR, regionDB = regionSetDB)
    # View(arrange(lResults, desc(oddsRatio)))
    # View(arrange(lResults, desc(pValueLog)))
    # sum(lResults$pValueLog < 2) / nrow(lResults)  
    
    return(lResults)
}

##################################
# running simulated data with noise
set.seed(1234)

noise = matrix(rnorm(n = nrow(both) * ncol(both), mean = 0, sd = 0.05), 
               nrow = nrow(both), ncol = ncol(both))
both = both +  noise
both[both < 0] = 0
both[both > 1] = 1
lResults = dmrLOLA(genomicSignal=both, signalCoord=signalCoord, 
                   GRList=GRList, rsAnno=rsAnno, targetVar=sampleStatus, 
                   dataID=paste0(dataID, "_gauss05"))
View(arrange(lResults, desc(oddsRatio)))
simpleCache(paste0("lolaResults", dataID, "_gauss05"), {
    lResults
}, assignToVariable = "lResults")

###################################
sampleN = 20
nCG = 10000
trueProp = 0.05
cgProps = rep(-1, nCG)
for (i in 1:nCG) {
  cgProps[i] = mean(sample(c(0, 1), size = sampleN, replace = TRUE, prob = c(1-trueProp, trueProp)))
}
hist(cgProps, breaks = seq(0, 1, by=0.05))
