# supervised COCOA for cancers other than BRCA


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

###############################################################################

# hist(pMeta$days_to_death)
# hist(pMeta$days_to_last_followup)
# pMeta = as.data.table(pMeta)
# a = pMeta[, .(mean(vital_status)), by=pathologic_stage]
# plot(c(2, 1, 3, 4, 0), a$V1)

trainDataInd = as.numeric(createDataPartition(allSampleLabels, p = (2/3), list=FALSE))

# give appropriate object names for downstream code
genomicSignal = methylMat[, trainDataInd]
sampleLabels = as.numeric(allSampleLabels[trainDataInd])
colsToAnnotate = "cancerStage"

# days = seq(365, 365* 10, by = 0.5)
# median(patientMetadata[patientMetadata$vital_status == "alive",]$days_to_last_follow_up, na.rm = TRUE)
# hist(patientMetadata[patientMetadata$vital_status == "dead",]$days_to_death)
# median(patientMetadata[patientMetadata$vital_status == "dead",]$days_to_death, na.rm = TRUE)
# brcaSurv = patientMetadata$vital_status
# sum(patientMetadata[patientMetadata$vital_status == "dead",]$days_to_death < 730)
# 
# 
# nAlive = nrow(patientMetadata[patientMetadata$vital_status == "alive",])
# totalDead = nrow(patientMetadata[patientMetadata$vital_status == "dead",])
# numberDead = (ecdf(x = patientMetadata[patientMetadata$vital_status == "dead",]$days_to_death)(days) * totalDead)
# numberKnownAlive = (1- ecdf(x = patientMetadata[patientMetadata$vital_status == "alive",]$days_to_last_follow_up)(days-1)) * nAlive
# plot(days, 
#      ecdf(x = patientMetadata[patientMetadata$vital_status == "dead",]$days_to_death )(days) * totalDead,  col = "red")
# lines(days, 
#       (1- ecdf(x = patientMetadata[patientMetadata$vital_status == "alive",]$days_to_last_follow_up)(days)) * nAlive, col="green")
# # ratio of those known to be alive to those known to be dead
# ratio = (numberKnownAlive + (totalDead - numberDead)) /
#     numberDead
# 
# plot(days, ratio)
# 
# daysCutoff = days[ratio <= 4][1]
# ratio[ratio <= 4][1]
# 

############################################################################
# run COCOA

simpleCache(paste0("rsScores_", dataID, "Cor"), {
    # create ATAC-protein correlation matrix
    actualCorMat = createCorFeatureMat(dataMat = genomicSignal,
                                       featureMat = as.matrix(sampleLabels),
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
# View(realRSScores[order(realRSScores$cancerStage, decreasing=TRUE), ])

##########################################################################
# permutation test for significance
# requires: nPerm, sampleLabels, genomicSignal, signalCoord, GRList, colsToAnnotate
# dataID
sampleLabels = data.frame(sampleLabels)
colnames(sampleLabels) = colsToAnnotate
source(ffProjCode("src/runPermTest.R"))

############################################################################

