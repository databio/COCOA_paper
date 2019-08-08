# permutation test by shuffling sample labels

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "19-applyPerm"
plotSubdir = "19-applyPerm/"

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 250

######################################################################
# required inputs to permutation test
dataID = "brcaDNAm657"

# assigns signalMat, signalCoord, loadingMat, pcScores
loadBRCADNAm()
genomicSignal = signalMat
sampleLabels = pcScores

# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg38")

colsToAnnotate = paste0("PC", 1:10)

variationMetric = "cov"

############################################################################

source(ffProjCode("src/runPermTest.R"))

############################################################################

load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", 
                                        dataID, ".RData")))
