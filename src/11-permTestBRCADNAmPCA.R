# permutation test by shuffling sample labels

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "11-permBRCA_DNAm_PCA"
plotSubdir = "11-permBRCA_DNAm_PCA/"

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 300
removeLowCov=TRUE
covCutoff = 100

######################################################################
# required inputs to permutation test
dataID = "brcaDNAm657"
colsToAnnotate = paste0("PC", 1:10)
variationMetric = "cov"

# assigns signalMat, signalCoord, loadingMat, pcScores
loadBRCADNAm()
genomicSignal = signalMat
sampleLabels = pcScores

# loads database of region sets 
# (assigns GRList, rsName, rsDescription, rsCollection to global environment)
loadGRList(genomeV="hg38")


############################################################################
# test median for reviewers
simpleCache(paste0("rsScores_", dataID, "_", variationMetric, "_median"), {
    #run COCOA
    actualResults = runCOCOA(genomicSignal =genomicSignal, 
                             signalCoord=signalCoord, GRList=GRList, 
                             signalCol = colsToAnnotate, 
                             variationMetric = variationMetric, 
                             targetVar = as.matrix(sampleLabels[, colsToAnnotate]),
                             scoringMetric = "regionMedian", verbose = TRUE)
    actualResults = cbind(actualResults, rsName=rsName, 
                          rsDescription=rsDescription, rsCollection=rsCollection)
    actualResults
}, assignToVariable = "realRSScores")

tmp = formattedCOCOAScores(rawScores = realRSScores, 
                           colsToAnnotate = colsToAnnotate, 
                           numTopRS = nrow(realRSScores))

write.csv(tmp, file = ffSheets(paste0("topRSScores","_", 
                                      dataID, "_", variationMetric, "_regionMedian.csv")), 
          row.names = FALSE)

############################################################################

simpleCache(paste0("rsScores_", dataID, "_", variationMetric), {
    #run COCOA
    actualResults = runCOCOA(genomicSignal =genomicSignal, 
                             signalCoord=signalCoord, GRList=GRList, 
                             signalCol = colsToAnnotate, 
                             variationMetric = variationMetric, 
                             targetVar = as.matrix(sampleLabels[, colsToAnnotate]),
                             scoringMetric = "regionMean", verbose = TRUE)
    actualResults = cbind(actualResults, rsName=rsName, 
                          rsDescription=rsDescription, rsCollection=rsCollection)
    actualResults
}, assignToVariable = "realRSScores")

tmp = formattedCOCOAScores(rawScores = realRSScores, 
                           colsToAnnotate = colsToAnnotate, 
                           numTopRS = nrow(realRSScores))

write.csv(tmp, file = ffSheets(paste0("topRSScores","_", 
                                      dataID, "_", variationMetric, ".csv")), 
          row.names = FALSE)

############################################################################

source(ffProjCode("runPermTest.R"))

############################################################################

# load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", nPerm, "_", variationMetric, 
#                    "_", dataID, ".RData")))
# load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", dataID, ".RData")))
conAssign(".analysisID", paste0("_", nPerm, "Perm_", variationMetric, "_", dataID))

# simpleCache(paste0("permPValsUncorrected", .analysisID), {
#     gPValDF
# }, recreate = FALSE, reload = TRUE, assignToVariable = "permPVals")
# 
# 
# formattedCOCOAScores(rawScores = realRSScores, colsToAnnotate = , 
#                      numTopRS = nrow(), pVals = , rankBy = )
