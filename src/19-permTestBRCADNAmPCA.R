# permutation test by shuffling sample labels

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

############################################################################

source(ffProjCode("src/runPermTest.R"))

############################################################################

load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", 
                                        dataID, ".RData")))
