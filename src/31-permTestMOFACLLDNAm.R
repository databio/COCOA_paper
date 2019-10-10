# permutation test by shuffling sample labels

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
overwriteRSScoreResultsCaches = FALSE

######################################################################
variationMetric = "cov"
colsToAnnotate = paste0("LF", 1:10)
# colsToAnnotate = paste0("LF", c(1:3, 5:7, 9))

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

#####################################################################
# permutation test for significance
# requires: nPerm, sampleLabels, genomicSignal, signalCoord, GRList, colsToAnnotate
# dataID, variationMetric (optional, default="cor")
source(ffProjCode("runPermTest.R"))

####################################################################

