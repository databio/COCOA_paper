# create survival plots

library(survival)
library(survminer)


# dataID = "kircMethyl"

sampleLabels = as.numeric(sampleLabels$cancerStage)

#############################################################################
simpleCache(paste0("rsScores_", dataID, "_", variationMetric), assignToVariable = "realRSScores")
simpleCache(paste0("pRankedScores", .analysisID), assignToVariable="pRankedScores")

topRSInd = rsRankingIndex(rsScores = pRankedScores, 
                          signalCol = list( "cancerStage_PValGroup", "cancerStage"), 
                          newColName="cancerStage")$cancerStage
orderedRSNames = rsName[topRSInd]
orderedRSDes = rsDescription[topRSInd]
thisRSInd = grep(pattern = "ezh2", x = orderedRSNames, ignore.case = TRUE)[1]
thisRSInd = c(thisRSInd, grep(pattern = "jund", x = orderedRSNames, ignore.case = TRUE)[1])
thisRSInd = c(thisRSInd, grep(pattern = "tcf7l2", x = orderedRSNames, ignore.case = TRUE)[1])
abbrevName = c("EZH2", "JUND", "TCF7L2")

actualCorMat = createCorFeatureMat(dataMat = genomicSignal,
                                   featureMat = as.matrix(sampleLabels),
                                   centerDataMat=TRUE, centerFeatureMat=TRUE,
                                   testType = variationMetric)
                                    
colnames(actualCorMat) <- colsToAnnotate

# look at average correlation per region to determine whether all CpGs 
# are going in same (important for aggregating since a simple average
# of the DNA methylation will run into problems if some CpGs go in opposite directions)
mByRegion = averagePerRegion(signal = actualCorMat, signalCoord = signalCoord, 
                             regionSet = GRList[[topRSInd[thisRSInd]]], 
                             signalCol = colsToAnnotate, absVal = FALSE)
hist(mByRegion$cancerStage)
# first test whether DNA methylation is associated with cancer stage 

#############################################################################
# "-" in name causes error
colnames(genomicSignal) = gsub(pattern = "-", replacement = "_", x = colnames(genomicSignal))
# average methylation per sample in these regions
mBySampleDF = runCOCOA(signal=genomicSignal, 
                    signalCoord=signalCoord, GRList=GRList[topRSInd[thisRSInd]], 
                    signalCol = colnames(genomicSignal), 
                    scoringMetric = "regionMean", verbose = TRUE)

# once for each region set
for (i in 1:nrow(mBySampleDF)) {
    mBySample = as.numeric(mBySampleDF[i, 1:ncol(genomicSignal)])
    trainMeta$methylScore = mBySample
    trainMeta$meanMethyl = colMeans(genomicSignal)
    hist(mBySample)
    plot(sampleLabels, mBySample)
    print(cor.test(x = sampleLabels, mBySample, method = "spearman"))
    print(cor.test(x = trainMeta$years_to_birth, mBySample))
    
    aliveMethyl = mBySample[trainMeta$vital_status == 0]
    deadMethyl = mBySample[trainMeta$vital_status == 1]
    print(wilcox.test(aliveMethyl, deadMethyl, conf.int = TRUE))
    
    # test whether DNA methylation is associated with survival
    ###### cox proportional hazards model 
    trainMeta$lastDate = trainMeta$days_to_last_followup
    trainMeta$lastDate[is.na(trainMeta$lastDate)] = trainMeta$days_to_death[is.na(trainMeta$lastDate)]
    
    # covariates
    covariateData = trainMeta
    patSurv = Surv(covariateData$lastDate, event=covariateData$vital_status)
    print(coxph(patSurv ~ years_to_birth + gender + meanMethyl + methylScore, data = covariateData))
    
    # kaplan meier
    # create groups
    covariateData$methylGroup = ecdf(covariateData$methylScore)(covariateData$methylScore)
    covariateData$methylGroup[covariateData$methylGroup < 0.25] = 0
    covariateData$methylGroup[covariateData$methylGroup > 0.75] = 2
    covariateData$methylGroup[(covariateData$methylGroup >= 0.25) & (covariateData$methylGroup <= 0.75)] = NA
    
    kmFit = survfit(patSurv ~ methylGroup, data=covariateData)
    ggsurvplot(kmFit)
}



################################################################################
# validation data
vGenomicSignal = methylMat[, -trainDataInd]
vSampleLabels = as.numeric(allSampleLabels[-trainDataInd])
colsToAnnotate = "cancerStage"

valMeta = pMeta[-trainDataInd, ]

# "-" in name causes error
colnames(vGenomicSignal) = gsub(pattern = "-", replacement = "_", x = colnames(vGenomicSignal))

mBySampleDF = runCOCOA(signal=vGenomicSignal, 
                     signalCoord=signalCoord, GRList=GRList[topRSInd[thisRSInd]], 
                     signalCol = colnames(vGenomicSignal), 
                     scoringMetric = "regionMean", verbose = TRUE)

for (i in 1:nrow(mBySampleDF)) {
    mBySample = as.numeric(mBySampleDF[i, 1:ncol(vGenomicSignal)])
    valMeta$methylScore = mBySample
    plot(vSampleLabels, mBySample)
    print(cor.test(x = vSampleLabels, mBySample, method="spearman"))
    print(cor.test(x = valMeta$years_to_birth, mBySample))
    
    aliveMethyl = mBySample[valMeta$vital_status == 0]
    deadMethyl = mBySample[valMeta$vital_status == 1]
    print(wilcox.test(aliveMethyl, deadMethyl, conf.int = TRUE))
    valMeta$meanMethyl = colMeans(vGenomicSignal)

    
    # test whether DNA methylation is associated with survival
    ###### cox proportional hazards model 
    valMeta$lastDate = valMeta$days_to_last_followup
    valMeta$lastDate[is.na(valMeta$lastDate)] = valMeta$days_to_death[is.na(valMeta$lastDate)]
    
    # covariates
    covariateData = valMeta
    patSurv = Surv(covariateData$lastDate, covariateData$vital_status)
    cModel = coxph(patSurv ~ years_to_birth + gender + meanMethyl + methylScore, data = covariateData)
    print(cModel)
    
    # create groups
    kmThresh = 0.25
    covariateData$methylGroup = ecdf(covariateData$methylScore)(covariateData$methylScore)
    covariateData$methylGroup[covariateData$methylGroup < kmThresh] = 0
    covariateData$methylGroup[covariateData$methylGroup > (1-kmThresh)] = 2
    covariateData$methylGroup[(covariateData$methylGroup >= kmThresh) & (covariateData$methylGroup <= (1-kmThresh))] = NA
    
    kmFit = survfit(patSurv ~ methylGroup, data=covariateData)
    kPlot = ggsurvplot(kmFit)
    kPlot
    pdf(ffPlot(paste0(plotSubdir, abbrevName[i], "_Kaplan", kmThresh, dataID, ".pdf")))
    print(kPlot)
    dev.off()
    # ggsave(filename = ffPlot(paste0(plotSubdir, "EZH2", "Kaplan", kmThresh, dataID, ".svg")) , plot = kPlot, device = "svg")
    
    plot(kmFit)
}


#############################################################################