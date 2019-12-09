# permutation test by shuffling sample labels

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

variationMetric = "cor"

# loads signalMat and signalCoord
loadBRCAatac()
genomicSignal = signalMat
simpleCache(paste0("brcaATACPCA_", ncol(genomicSignal)))
# the PC scores
sampleLabels = brcaATACPCA_73$x

dataID = paste0("brcaATAC", ncol(genomicSignal))

colsToAnnotate = paste0("PC", 1:10)
AGR2
C10orf10
CDS1
FAM169A
FXYD3
KLK10
LAD1
MTUS1
PLS1
PRRG4
RHOD
SERPINB1
SLPI
TMEM30B
TPD52L1
TSPAN1
ZHX2
ZNF165

### get shared samples and put data in same order 
sharedSamples = colnames(signalMat)[colnames(signalMat) %in% row.names(sampleLabels)]
genomicSignal = signalMat[, sharedSamples]
sampleLabels = sampleLabels[sharedSamples, colsToAnnotate]

# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg38")


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
}, assignToVariable = "realRSScores", recreate=TRUE)

############################################################################

# requires: nPerm, sampleLabels, genomicSignal, signalCoord, GRList, colsToAnnotate
# dataID
source(ffProjCode("runPermTest.R"))

############################################################################

# load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", nPerm, "Perm_", variationMetric, 
#                    "_", dataID, ".RData")))
# rsPermScores = ret

# function requirements:
# loads existing object
# starts script/code at specified point
# adds new object to existing object
# resume=TRUE simpleCache
# 