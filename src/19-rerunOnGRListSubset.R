# running COCOA perm test on new region sets and adding to existing results

###############################################################################
# functions to append new caches to old caches


# append realRSscores
appendRealRSScores = function(realRSScores1, realRSScores2) {
    if (ncol(realRSScores1) != ncol(realRSScores2)) {
        stop()
    }
    # assuming more recent col names are correct
    if (any(colnames(realRSScores1) != colnames(realRSScores2))) {
        colnames(realRSScores1) = colnames(realRSScores2)
    }
    
    return(rbind(realRSScores1, realRSScores2))
}

# append rsPermScores
# @param rsPermScores1 list.
appendRSPermScores = function(rsPermScores1, rsPermScores2) {
    if (ncol(rsPermScores1[[1]]) != ncol(rsPermScores2[[1]])) {
        stop()
    }
    # assuming more recent col names are correct
    if (any(colnames(rsPermScores1[[1]]) != colnames(rsPermScores2[[1]]))) {
        for (i in seq_along(rsPermScores1)) {
            colnames(rsPermScores1[[i]]) = colnames(rsPermScores2[[i]])
        }
    }
    
    
    combRSPermScores = list()
    # add to each list item
    for (i in seq_along(rsPermScores1)) {
        combRSPermScores[[i]] = rbind(rsPermScores1[[i]], rsPermScores2[[i]])
    }
    return(combRSPermScores)
}


#############################################################################
# BRCA DNA methylation PCA
#############################################################################
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "19-permBRCA_DNAm_PCA"
plotSubdir = "19-permBRCA_DNAm_PCA/"

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 300

######################################################################
# required inputs to permutation test
dataID = "brcaDNAm657"
devtools::load_all(ffCode("COCOA/"))
colsToAnnotate = paste0("PC", 1:10)
variationMetric = "cov"

# assigns signalMat, signalCoord, loadingMat, pcScores
loadBRCADNAm()
genomicSignal = signalMat
sampleLabels = pcScores

# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg38")

# just run on this subset
dataID = paste0(dataID, "_hemaATAC")
hemaATACInd = rsCollection == "hematopoietic_ATACseq_GSE75384"
GRList = GRList[hemaATACInd]
rsName = rsName[hemaATACInd]
rsDescription = rsDescription[hemaATACInd]
rsCollection = rsCollection[hemaATACInd]

############################################################################


                
simpleCache(paste0("rsScores_", dataID, "_", variationMetric), {
    # create ATAC-protein correlation matrix
    actualCorMat = createCorFeatureMat(dataMat = genomicSignal,
                                       featureMat = as.matrix(sampleLabels[, colsToAnnotate]),
                                       centerDataMat=TRUE, centerFeatureMat=TRUE, testType = variationMetric)
    colnames(actualCorMat) <- colsToAnnotate
    
    #run COCOA
    actualResults = runCOCOA(signal=actualCorMat, 
                             signalCoord=signalCoord, GRList=GRList, 
                             signalCol = colsToAnnotate, 
                             scoringMetric = "default", verbose = TRUE)
    actualResults = cbind(actualResults, rsName=rsName, 
                          rsDescription=rsDescription, rsCollection=rsCollection)
    actualResults
}, assignToVariable = "realRSScores")

simpleCache(paste0("rsScores_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID), 
                   "_", variationMetric), 
            assignToVariable = "realRSScores1")

simpleCache(paste0("rsScores_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID), 
                   "_", variationMetric), {    
                       combRealRSScores = appendRealRSScores(realRSScores1 = realRSScores1, realRSScores2 = realRSScores)
                       combRealRSScores
        }, recreate=TRUE)



############################################################################

source(ffProjCode("runPermTest.R"))

simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), assignToVariable = "rsPermScores")
simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID)), assignToVariable = "rsPermScores1")

simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID)), {
                       combRSPermScores = appendRSPermScores(rsPermScores1 = rsPermScores1, rsPermScores2 = rsPermScores)
                       combRSPermScores
                   }, recreate = TRUE)

############################################################################

# load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", nPerm, "_", variationMetric, 
#                    "_", dataID, ".RData")))
# load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", dataID, ".RData")))

#############################################################################
# BRCA ATAC PCA
#############################################################################

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
devtools::load_all(ffCode("COCOA"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "19-permTestBRCA_ATAC_PCA"
plotSubdir = "19-permBRCA_ATAC_PCA/"
sheetsDir = ffProc("COCOA_paper/analysis/sheets/")

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 300

######################################################################
# required inputs to permutation test
dataID = "brcaATACPCAPerm"
variationMetric = "cor"

# loads signalMat and signalCoord
loadBRCAatac()
genomicSignal = signalMat
simpleCache(paste0("brcaATACPCA_73"))
# the PC scores
sampleLabels = brcaATACPCA_73$x

colsToAnnotate = paste0("PC", 1:10)

### get shared samples and put data in same order 
sharedSamples = colnames(signalMat)[colnames(signalMat) %in% row.names(sampleLabels)]
genomicSignal = signalMat[, sharedSamples]
sampleLabels = sampleLabels[sharedSamples, colsToAnnotate]

# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg38")


# just run on this subset
dataID = paste0(dataID, "_hemaATAC")
hemaATACInd = rsCollection == "hematopoietic_ATACseq_GSE75384"
GRList = GRList[hemaATACInd]
rsName = rsName[hemaATACInd]
rsDescription = rsDescription[hemaATACInd]
rsCollection = rsCollection[hemaATACInd]

######################################################################

# get the "true" COCOA scores before doing the permutation test
simpleCache(paste0("rsScores_", dataID, "_", variationMetric), {
    # create ATAC-protein correlation matrix
    actualCorMat = createCorFeatureMat(dataMat = genomicSignal,
                                       featureMat = as.matrix(sampleLabels[, colsToAnnotate]),
                                       centerDataMat=TRUE, centerFeatureMat=TRUE)
    colnames(actualCorMat) <- colsToAnnotate
    
    #run COCOA
    actualResults = runCOCOA(signal=actualCorMat, 
                             signalCoord=signalCoord, GRList=GRList, 
                             signalCol = colsToAnnotate, 
                             scoringMetric = "default", verbose = TRUE)
    actualResults = cbind(actualResults, rsName=rsName, 
                          rsDescription=rsDescription)
    actualResults
}, assignToVariable = "realRSScores")

simpleCache(paste0("rsScores_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID), 
                   "_", variationMetric), 
            assignToVariable = "realRSScores1")

simpleCache(paste0("rsScores_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID), 
                   "_", variationMetric), {    
                       combRealRSScores = appendRealRSScores(realRSScores1 = realRSScores1, realRSScores2 = realRSScores)
                       combRealRSScores
                   }, recreate=TRUE)

############################################################################

# requires: nPerm, sampleLabels, genomicSignal, signalCoord, GRList, colsToAnnotate
# dataID
source(ffProjCode("runPermTest.R"))

simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), assignToVariable = "rsPermScores")
simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID)), assignToVariable = "rsPermScores1")

simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID)), {
                       combRSPermScores = appendRSPermScores(rsPermScores1 = rsPermScores1, rsPermScores2 = rsPermScores)
                       combRSPermScores
                   }, recreate = TRUE)


#############################################################################
# MOFA DNA methylation
#############################################################################


source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
devtools::load_all(ffCode("COCOA/"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "19-permMOFACLLDNAm"
plotSubdir = "19-permMOFACLLDNAm/"
dataID = "CLL196MOFA"
sheetsDir = ffProc("COCOA_paper/analysis/sheets/")

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 300

######################################################################
variationMetric = "cov"
# simpleCache(paste0("rsScore_", paste0(dataID, "_", variationMetric)), assignToVariable = "realRSScores")
colsToAnnotate = paste0("LF", 1:10)

# assigns methylMat, signalCoord, latentFactorMat to global environment
loadMOFAData()
genomicSignal = methylMat
sampleLabels = data.frame(latentFactorMat)

# loads database of region sets 
# (assigns GRList, rsName, rsDescription, rsCollection to global environment)
loadGRList(genomeV="hg19")

# row.names(realRSScores) = realRSScores$rsName
# sharedRSNames = names(GRList)[names(GRList) %in% realRSScores$rsName]
# GRList = GRList[sharedRSNames]
# realRSScores = realRSScores[sharedRSNames, ]

# colsToAnnotate = paste0("LF", c(1:3, 5:7, 9))

# just run on this subset
dataID = paste0(dataID, "_hemaATAC")
hemaATACInd = rsCollection == "hematopoietic_ATACseq_GSE75384"
GRList = GRList[hemaATACInd]
rsName = rsName[hemaATACInd]
rsDescription = rsDescription[hemaATACInd]
rsCollection = rsCollection[hemaATACInd]

rsScoreCacheName = paste0("rsScore_", paste0(dataID, "_", variationMetric))

############################################################################
# create the actual MOFA scores

simpleCache(rsScoreCacheName, {
    
    # calculate correlation/covariation with each latent factor from MOFA
    featurePCCor = createCorFeatureMat(dataMat = genomicSignal, 
                                       featureMat = sampleLabels, 
                                       centerDataMat = TRUE, 
                                       centerFeatureMat = TRUE, 
                                       testType = variationMetric)
    # run COCOA analysis
    
    rsScore = runCOCOA(signal=abs(featurePCCor), 
                       signalCoord = signalCoord, 
                       GRList, 
                       signalCol = colsToAnnotate, 
                       signalCoordType = "singleBase")
    rsScore$rsName = rsName
    rsScore$rsDescription= rsDescription
    rsScore$rsCollection = rsCollection
    rsScore
}, recreate=overwriteRSScoreResultsCaches, assignToVariable = "realRSScores")
# rsScores = as.data.table(rsScore)

simpleCache(paste0("rsScores_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID), 
                   "_", variationMetric), 
            assignToVariable = "realRSScores1")

simpleCache(paste0("rsScores_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID), 
                   "_", variationMetric), {    
                       combRealRSScores = appendRealRSScores(realRSScores1 = realRSScores1, realRSScores2 = realRSScores)
                       combRealRSScores
                   }, recreate=TRUE)


#####################################################################
# permutation test for significance
# requires: nPerm, sampleLabels, genomicSignal, signalCoord, GRList, colsToAnnotate
# dataID, variationMetric (optional, default="cor")
source(ffProjCode("runPermTest.R"))


#############################################################################
# KIRC methylation
#############################################################################


source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
devtools::load_all(ffCode("COCOA"))
library(caret)

scriptID = "22-tcgaCorCOCOA"
plotSubdir = "22-tcgaCorCOCOA/"
sheetsDir = ffProc("COCOA_paper/analysis/sheets/")
dataID = "kircMethyl"

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 300

variationMetric = "spearmanCor" 
######################################################################
# load data

# loads methylList, pMeta (patient metadata)
loadTCGAMethylation(cancerID = "KIRC")
patientMetadata = pMeta
methylMat = methylList$methylProp
signalCoord = methylList$coordinates

sampleType = substr(colnames(methylMat), start = 14, stop = 15)
# 01 is primary solid tumor, 11 is solid normal tissue, 05 is new primary tumor
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
normalSampleInd = (sampleType == "11")
tumorSampleInd = (sampleType == "01") # exclude the extra sample 05
methylMat = methylMat[, tumorSampleInd]
# now I only need patient ID
colnames(methylMat) = substr(colnames(methylMat), start = 1, stop = 12)

## order samples consistently
pMeta = pMeta[colnames(methylMat), ]
# screen out patients without stage
naInd = is.na(pMeta$pathologic_stage)
methylMat = methylMat[, !naInd]
pMeta = pMeta[!naInd, ]

allSampleLabels = factor(pMeta$pathologic_stage, levels = c("stage i", "stage ii", "stage iii", "stage iv"))

loadGRList(genomeV = "hg19")

# just run on this subset
dataID = paste0(dataID, "_hemaATAC")
hemaATACInd = rsCollection == "hematopoietic_ATACseq_GSE75384"
GRList = GRList[hemaATACInd]
rsName = rsName[hemaATACInd]
rsDescription = rsDescription[hemaATACInd]
rsCollection = rsCollection[hemaATACInd]

###############################################################################

trainDataInd = as.numeric(createDataPartition(allSampleLabels, p = (2/3), list=FALSE))

# give appropriate object names for downstream code
genomicSignal = methylMat[, trainDataInd]
sampleLabels = as.numeric(allSampleLabels[trainDataInd])
colsToAnnotate = "cancerStage"

trainMeta = pMeta[trainDataInd, ]

############################################################################
# run COCOA

simpleCache(paste0("rsScores_", dataID, "_", variationMetric), {
    # create ATAC-protein correlation matrix
    actualCorMat = createCorFeatureMat(dataMat = genomicSignal,
                                       featureMat = as.matrix(sampleLabels),
                                       centerDataMat=TRUE, centerFeatureMat=TRUE, testType = variationMetric)
    colnames(actualCorMat) <- colsToAnnotate
    
    #run COCOA
    actualResults = runCOCOA(signal=actualCorMat, 
                             signalCoord=signalCoord, GRList=GRList, 
                             signalCol = colsToAnnotate, 
                             scoringMetric = "default", verbose = TRUE)
    actualResults = cbind(actualResults, rsName=rsName, 
                          rsDescription=rsDescription)
    actualResults
}, assignToVariable = "realRSScores")
# View(realRSScores[order(realRSScores$cancerStage, decreasing=TRUE), ])

simpleCache(paste0("rsScores_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID), 
                   "_", variationMetric), 
            assignToVariable = "realRSScores1")

simpleCache(paste0("rsScores_", 
                   sub(pattern = "_hemaATAC", replacement = "", x = dataID), 
                   "_", variationMetric), {    
                       combRealRSScores = appendRealRSScores(realRSScores1 = realRSScores1, realRSScores2 = realRSScores)
                       combRealRSScores
                   }, recreate=TRUE)

##########################################################################
# permutation test for significance
# requires: nPerm, sampleLabels, genomicSignal, signalCoord, GRList, colsToAnnotate
# dataID
sampleLabels = data.frame(sampleLabels)
colnames(sampleLabels) = colsToAnnotate

source(ffProjCode("runPermTest.R"))

###############################################################################




