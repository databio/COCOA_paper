# permutation test by shuffling sample labels

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
devtools::load_all(ffCode("COCOA"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "19-permTestBRCA_ATAC_Protein"
plotSubdir = "19-permBRCA_ATAC_Protein/"
sheetsDir = ffProc("COCOA_paper/analysis/sheets/")

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 250

######################################################################
# load protein data

mainP = c("NFKBP65_pS536", "GATA3", "GATA6", "ERALPHA_pS118", "ERALPHA", "AR", 
          "CJUN_pS73", "CMYC", "BETACATENIN", "PR", "SMAD1", "SMAD3", "STAT3_pY705", 
          "STAT5ALPHA", "ETS1", "FOXM1", "IRF1")
alsoOfInterest = c("CD31", "BRAF", "CYCLIND1", "RICTOR", "CD49B", "CD20")
proteinOfInterest = c(mainP, alsoOfInterest)

prot = read.csv(file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/TCGA-BRCA-L4.csv"))
prot$subject_ID = substr(prot$Sample_ID, start = 1, stop = 12)

# figuring out which patients are represented more than once
idt = table(prot$subject_ID)
# head(sort(idt, decreasing = TRUE))
multiInd = which(idt > 1)
multiPatient = names(multiInd)

# remove patients with more than one sample
prot = prot[!(prot$subject_ID %in% multiPatient), ]
row.names(prot) = prot$subject_ID

######################################################################
# required inputs to permutation test



# loads signalMat and signalCoord
loadBRCAatac()

### get shared samples and put data in same order 
sharedSamples = colnames(signalMat)[colnames(signalMat) %in% prot$subject_ID]
genomicSignal = signalMat[, sharedSamples]
sampleLabels = prot[sharedSamples, proteinOfInterest]

# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg38")

colsToAnnotate = proteinOfInterest #protein columns

dataID = "brcaATACProteinPerm"


# get the "true" COCOA scores before doing the permutation test
simpleCache("rsScores_brcaProteinATACCor", {
    # create ATAC-protein correlation matrix
    actualCorMat = createCorFeatureMat(dataMat = genomicSignal,
                                       featureMat = as.matrix(sampleLabels[, proteinOfInterest]),
                                       centerDataMat=TRUE, centerFeatureMat=TRUE)
    colnames(actualCorMat) <- proteinOfInterest
    
    #run COCOA
    actualResults = runCOCOA(signal=actualCorMat, 
                           signalCoord=signalCoord, GRList=GRList, 
                           signalCol = colsToAnnotate, 
                           scoringMetric = "default", verbose = TRUE)
    actualResults = cbind(actualResults, rsName=rsName, 
                          rsDescription=rsDescription)
    actualResults
}, assignToVariable = "realRSScores")





############################################################################

# 
# multiNiceHist(file = ffPlot(paste0(plotSubdir,"proteinATACCorDist.pdf")), dataDF = as.data.frame(trueCorMat),
#               colsToPlot = proteinOfInterest, xLabels = "Correlation of protein with ATAC-seq counts", binwidth=0.05, boundary=0,
#               plotTitles = proteinOfInterest, yLabel = "Number of regions")

# requires: nPerm, sampleLabels, genomicSignal, signalCoord, GRList, colsToAnnotate
# dataID
source(ffProjCode("src/runPermTest.R"))

############################################################################




