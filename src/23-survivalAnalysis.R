# create survival plots

library(survival)
library(survminer)


dataID = "kircMethyl"

#############################################################################
simpleCache(paste0("rsScores_", dataID, "Cor"), assignToVariable = "realRSScores")

topRSInd = rsRankingIndex(rsScores = realRSScores, signalCol = "cancerStage")$cancerStage

actualCorMat = createCorFeatureMat(dataMat = genomicSignal,
                                   featureMat = as.matrix(sampleLabels),
                                   centerDataMat=TRUE, centerFeatureMat=TRUE)
colnames(actualCorMat) <- colsToAnnotate

# look at average correlation per region to determine whether all CpGs 
# are going in same (important for aggregating since a simple average
# of the DNA methylation will run into problems if some CpGs go in opposite directions)
mByRegion = averagePerRegion(signal = actualCorMat, signalCoord = signalCoord, 
                             regionSet = GRList[[topRSInd[1]]], 
                             signalCol = colsToAnnotate, absVal = FALSE)
hist(mByRegion$cancerStage)
# first test whether DNA methylation is associated with cancer stage 

#############################################################################
# "-" in name causes error
colnames(genomicSignal) = gsub(pattern = "-", replacement = "_", x = colnames(genomicSignal))
# average methylation per sample in these regions
mBySample = runCOCOA(signal=genomicSignal, 
                    signalCoord=signalCoord, GRList=GRList[topRSInd[1]], 
                    signalCol = colnames(genomicSignal), 
                    scoringMetric = "regionMean", verbose = TRUE)
mBySample = as.numeric(mBySample[1:ncol(genomicSignal)])
trainMeta$methylScore = mBySample
trainMeta$meanMethyl = colMeans(genomicSignal)
hist(mBySample)
plot(sampleLabels, mBySample)
cor.test(x = sampleLabels, mBySample)
cor.test(x = trainMeta$years_to_birth, mBySample)

aliveMethyl = mBySample[trainMeta$vital_status == 0]
deadMethyl = mBySample[trainMeta$vital_status == 1]
wilcox.test(aliveMethyl, deadMethyl)
wilcox.test(aliveMethyl, deadMethyl, conf.int = TRUE)


################################################################################
# validation data
vGenomicSignal = methylMat[, -trainDataInd]
vSampleLabels = as.numeric(allSampleLabels[-trainDataInd])
colsToAnnotate = "cancerStage"

valMeta = pMeta[-trainDataInd, ]

# "-" in name causes error
colnames(vGenomicSignal) = gsub(pattern = "-", replacement = "_", x = colnames(vGenomicSignal))

mBySample = runCOCOA(signal=vGenomicSignal, 
                     signalCoord=signalCoord, GRList=GRList[topRSInd[1]], 
                     signalCol = colnames(vGenomicSignal), 
                     scoringMetric = "regionMean", verbose = TRUE)

mBySample = as.numeric(mBySample[1:ncol(vGenomicSignal)])
valMeta$methylScore = mBySample
plot(vSampleLabels, mBySample)
cor.test(x = vSampleLabels, mBySample)
cor.test(x = valMeta$years_to_birth, mBySample)

aliveMethyl = mBySample[valMeta$vital_status == 0]
deadMethyl = mBySample[valMeta$vital_status == 1]
wilcox.test(aliveMethyl, deadMethyl, conf.int = TRUE)
valMeta$meanMethyl = colMeans(vGenomicSignal)

################################################################################
# test whether DNA methylation is associated with survival
###### cox proportional hazards model 
trainMeta$lastDate = trainMeta$days_to_last_followup
trainMeta$lastDate[is.na(trainMeta$lastDate)] = trainMeta$days_to_death[is.na(trainMeta$lastDate)]

# covariates
covariateData = trainMeta
patSurv = Surv(covariateData$lastDate, event=covariateData$vital_status)
coxph(patSurv ~ years_to_birth + gender + meanMethyl + methylScore, data = covariateData)

# kaplan meier
# create groups
covariateData$methylGroup = ecdf(covariateData$methylScore)(covariateData$methylScore)
covariateData$methylGroup[covariateData$methylGroup < 0.25] = 0
covariateData$methylGroup[covariateData$methylGroup > 0.75] = 2
covariateData$methylGroup[(covariateData$methylGroup >= 0.25) & (covariateData$methylGroup <= 0.75)] = NA

kmFit = survfit(patSurv ~ methylGroup, data=covariateData)
ggsurvplot(kmFit)
plot(kmFit)

###
valMeta$lastDate = valMeta$days_to_last_followup
valMeta$lastDate[is.na(valMeta$lastDate)] = valMeta$days_to_death[is.na(valMeta$lastDate)]

# covariates
covariateData = valMeta
patSurv = Surv(covariateData$lastDate, covariateData$vital_status)
cModel = coxph(patSurv ~ years_to_birth + gender + meanMethyl + methylScore, data = covariateData)

# create groups
kmThresh = 0.25
covariateData$methylGroup = ecdf(covariateData$methylScore)(covariateData$methylScore)
covariateData$methylGroup[covariateData$methylGroup < kmThresh] = 0
covariateData$methylGroup[covariateData$methylGroup > (1-kmThresh)] = 2
covariateData$methylGroup[(covariateData$methylGroup >= kmThresh) & (covariateData$methylGroup <= (1-kmThresh))] = NA

kmFit = survfit(patSurv ~ methylGroup, data=covariateData)
kPlot = ggsurvplot(kmFit)
svg(ffPlot(paste0(plotSubdir, "EZH2", "Kaplan", kmThresh, dataID, ".svg")))
kPlot
dev.off()
# ggsave(filename = ffPlot(paste0(plotSubdir, "EZH2", "Kaplan", kmThresh, dataID, ".svg")) , plot = kPlot, device = "svg")

plot(kmFit)


#########

# kaplan meier plot: groups are predicted good outcome vs predicted bad outcome
patSurv = Surv(patientMetadata_pqc$os_months, rep(1, length(patientMetadata_pqc$os_months)))
kmFit = survfit(patSurv ~ 1)
patSurv2 = Surv(patientMetadata_pqc$os_months, rep_len(c(1,1,0), length.out = length(patientMetadata_pqc$os_months)))
kmFit2 = survfit(formula = patSurv2 ~ patientMetadata_pqc$Complex, data = patientMetadata_pqc)

plot(kmFit2)


