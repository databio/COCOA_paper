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

trainMeta = pMeta[trainDataInd, ]

############################################################################
# # test whether cancer stages have genomewide differences in DNA methylation levels
# sampleMeanMethyl = colMeans(genomicSignal)
# sampleMeanDT = data.table(cancerStage= sampleLabels, meanMethyl=sampleMeanMethyl)
# meanByStage = sampleMeanDT[, .(meanStageMethyl = mean(meanMethyl)), by=cancerStage]
# plot(sampleMeanDT$cancerStage, sampleMeanDT$meanMethyl)
# cor.test(sampleMeanDT$cancerStage, sampleMeanDT$meanMethyl)
# 
# # normalize for average methylation level, by cancer stage
# for (i in seq_along(unique(sampleLabels))) {
#     
#     # normalize one cancer stage
#     genomicSignal[, sampleLabels == i] = genomicSignal[, sampleLabels == i] - meanByStage$meanStageMethyl[i]
#     
# }
# 
# partCorMat = createCorFeatureMat(dataMat = genomicSignal,
#                                    featureMat = as.matrix(sampleLabels),
#                                    centerDataMat=FALSE, centerFeatureMat=TRUE, testType = "pcor",
#                                  covariate = sampleMeanMethyl)
# colnames(partCorMat) <- colsToAnnotate
# 
# #run COCOA
# partCorResults = runCOCOA(signal=partCorMat, 
#                          signalCoord=signalCoord, GRList=GRList, 
#                          signalCol = colsToAnnotate, 
#                          scoringMetric = "default", verbose = TRUE)
# partCorResults = cbind(partCorResults, rsName=rsName, 
#                       rsDescription=rsDescription)
# 
# dataID = paste0(dataID, "Norm")
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

##########################################################################
# permutation test for significance
# requires: nPerm, sampleLabels, genomicSignal, signalCoord, GRList, colsToAnnotate
# dataID
sampleLabels = data.frame(sampleLabels)
colnames(sampleLabels) = colsToAnnotate

source(ffProjCode("runPermTest.R"))


