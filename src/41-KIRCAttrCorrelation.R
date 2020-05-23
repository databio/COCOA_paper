# supervised COCOA for cancers other than BRCA


source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(COCOA)
library(caret)
library(minfi)

scriptID = "41-tcgaCorCOCOA"
plotSubdir = "41-tcgaCorCOCOA/"
sheetsDir = ffProc("COCOA_paper/analysis/sheets/")


if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 300
dataID = "kircMethyl"
variationMetric = "spearmanCor" 
removeLowCov = TRUE
covCutoff = 100

######################################################################
# load data

# loads the following variables to current environment:
# c("methylMat", "signalCoord", "pMeta", "allSampleLabels"),
loadProcessKIRCMethyl()

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
sampleLabels = data.frame(sampleLabels)
colnames(sampleLabels) = colsToAnnotate

trainMeta = pMeta[trainDataInd, ]
dataID = paste0(dataID, ncol(genomicSignal))

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
# dataID, optional: recreate (for simpleCache)

source(ffProjCode("runPermTest.R"))

##########################################################################
# differential methylation analysis of cancer stage
# for comparison with COCOA

# separate into two groups
# stage1M = genomicSignal[, sampleLabels <= 1]
# stage234M = genomicSignal[, sampleLabels >= 2]

# try modeling as continuous variable
designMat = model.matrix(object = ~ cancerStage, data=sampleLabels)
# formula(~ cancerStage, data=sampleLabels)

# diff methyl CpGs
dmp <- dmpFinder(genomicSignal, pheno = sampleLabels$cancerStage, type = "continuous")
hist(-log10(dmp$qval))

# find differentially methylated regions
signalCoord$index = 1:nrow(signalCoord)
signalCoord = arrange(signalCoord, chr, start)
genomicSignal = genomicSignal[signalCoord$index, ]
# figure out how CpGs are being grouped and if altered CpG order is meesing things up
a =bumphunter(object = genomicSignal, design=designMat, chr=signalCoord$chr, 
           pos=signalCoord$start, coef=2, cutoff = 0.01, B=10)# , type="Beta")
any(table(a$table$cluster) > 1)
table(a$table$cluster)[table(a$table$cluster) > 1]
a$table[names(table(a$table$cluster)[table(a$table$cluster) > 1]), ]

# use LOLA to identify region sets associated with
LOLA::runLOLA()


