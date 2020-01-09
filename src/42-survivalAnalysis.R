# create survival plots
# ideas for figures and what information is included in survival plots:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6245494/
# also has table of number of patients at risk, no tick marks on line but has confidence interval

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3008568/
# didn't inlcude on plot how many patients were in each group but had censorship marks on lines
# has table in figure for Cox PH model, columns: Variables in Model, p value, Adjusted Hazard Ratio
# 95 % CI for HR (Lower and Upper cols)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6676410/


library(survival)
library(survminer)
library(tidyr)
plotWidth = 100
plotHeight = 100
plotUnit = "mm"
if (!exists("dataID")) {
    dataID = "kircMethyl214"
}

if (!exists("variationMetric")) {
    variationMetric = "spearmanCor"    
}
if (!exists("nPerm")) {
    nPerm = 300
}
if (!exists(".analysisID")) {
    .analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)
}

################################################################################

simpleCache(paste0("rsScores_", dataID, "_", variationMetric), assignToVariable = "realRSScores")
keepInd = realRSScores$regionSetCoverage >= 100
realRSScores = realRSScores[keepInd, ]
GRList = GRList[keepInd]
rsName = rsName[keepInd]
rsDescription = rsDescription[keepInd]
rsCollection = rsCollection[keepInd]

###############################################################################
sampleLabels = as.numeric(sampleLabels$cancerStage)

#############################################################################

# simpleCache(paste0("pRankedScores", .analysisID), assignToVariable="pRankedScores")

# topRSInd = rsRankingIndex(rsScores = pRankedScores, 
#                           signalCol = list( "cancerStage_PValGroup", "cancerStage"), 
#                           newColName="cancerStage")$cancerStage
# orderedRSNames = pRankedScores$rsName[topRSInd]
# orderedRSDes = pRankedScores$rsDescription[topRSInd]
topRSInd = rsRankingIndex(rsScores = realRSScores,
                          signalCol = "cancerStage",
                          newColName="cancerStage")$cancerStage
orderedRSNames = realRSScores$rsName[topRSInd]
orderedRSDes = realRSScores$rsDescription[topRSInd]


# names to search for 
topRSSearchNames = c("ezh2", "jund", "tcf7l2")
# topRSSearchNames = c("ezh2", "jund", "tcf7l2", "stemNotDiff", "P300", "Ctbp2", "H1hesc_E")
thisRSInd = grep(pattern = "ezh2", x = orderedRSNames, ignore.case = TRUE)[1]
for (i in 2:length(topRSSearchNames)) {
    thisRSInd = c(thisRSInd, grep(pattern = topRSSearchNames[i], x = orderedRSNames, ignore.case = TRUE)[1])
}

abbrevName = c("EZH2", "JUND", "TCF7L2")
# abbrevName = c("EZH2", "JUND", "TCF7L2", "stemNotDiff", "P300", "Ctbp2", "H1hesc_E")
# interesting: stemNotDiff, P300, Ctbp2, H1hesc_E

# actualCorMat = COCOA:::createCorFeatureMat(dataMat = genomicSignal,
#                                    featureMat = as.matrix(sampleLabels),
#                                    centerDataMat=TRUE, centerFeatureMat=TRUE,
#                                    testType = variationMetric)
#                                     
# colnames(actualCorMat) <- colsToAnnotate
# 
# # look at average correlation per region to determine whether all CpGs 
# # are going in same (important for aggregating since a simple average
# # of the DNA methylation will run into problems if some CpGs go in opposite directions)
# mByRegion = COCOA:::averagePerRegion(signal = actualCorMat, signalCoord = signalCoord, 
#                              regionSet = GRList[[topRSInd[1]]], 
#                              signalCol = colsToAnnotate, absVal = FALSE)
# hist(mByRegion$cancerStage)
# sum(mByRegion$cancerStage < 0) / length(mByRegion$cancerStage)

#############################################################################
# KIRC Panel C
# first test whether DNA methylation is associated with cancer stage 


# "-" in name causes error
colnames(genomicSignal) = gsub(pattern = "-", replacement = "_", x = colnames(genomicSignal))
# average methylation per sample in these regions
mBySampleDF = aggregateSignalGRList(signal=genomicSignal, 
                    signalCoord=signalCoord, GRList=GRList[topRSInd[thisRSInd]], 
                    signalCol = colnames(genomicSignal), 
                    scoringMetric = "regionMean", verbose = TRUE)
# convert to long format for plotting
longMByS = transpose(mBySampleDF) 
colnames(longMByS) = abbrevName
longMByS$subjectID = colnames(mBySampleDF)
longMByS = longMByS[1:ncol(genomicSignal), ]
longMByS = cbind(longMByS, cancerStage=as.factor(sampleLabels))
longMByS = gather(data = longMByS, "regionSetName", "methylScore", abbrevName)
mByStagePlot = ggplot(data=longMByS, mapping = aes(x=cancerStage, y=methylScore)) + geom_violin() + facet_wrap("regionSetName") + 
    geom_smooth(mapping = aes(x = as.numeric(cancerStage), y=methylScore)) + xlab("Cancer stage") + ylab("Methylation proportion")
mByStagePlot
ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStageTraining.svg")), 
       plot = mByStagePlot, device = "svg")

# individual RS plots
for (i in seq_along(abbrevName)) {
    mByStagePlotSingle = ggplot(data=filter(longMByS, regionSetName == abbrevName[i]), 
                              mapping = aes(x=cancerStage, y=methylScore)) + 
        geom_violin() +
        geom_smooth(mapping = aes(x = as.numeric(cancerStage), y=methylScore)) + xlab("Cancer stage") + ylab("Methylation proportion")
    mByStagePlotSingle
    ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStageTraining", abbrevName[i], ".svg")), 
           plot = mByStagePlotSingle, device = "svg", width=plotWidth, height = plotHeight / 2,
           units = plotUnit)
}




# test how much methylation level in each region set is correlated with
# methylation level in other region sets
corMat = cor(t(mBySampleDF[, 1:ncol(genomicSignal)]))

###############################################################################
# Panel

# data.frame to store results 
tmp = rep(-999, nrow(mBySampleDF))
trainResDF = data.frame(corPVal=tmp, 
                        corCoef=tmp, 
                        coxPVal=tmp, 
                        coxEffect=tmp, 
                        coxModel=tmp)

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
    patSurv = Surv(covariateData$lastDate / (365/12), event=covariateData$vital_status)
    myModel = coxph(patSurv ~ years_to_birth + gender + meanMethyl + methylScore, data = covariateData)
    
    sink(file = ffPlot(paste0(plotSubdir, "coxphModel", abbrevName[i], ".txt")))
        print(myModel)
        print(summary(myModel))
    sink()
    
        
    # kaplan meier
    # create groups
    covariateData$methylGroup = ecdf(covariateData$methylScore)(covariateData$methylScore)
    covariateData$methylGroup[covariateData$methylGroup < 0.25] = 0
    covariateData$methylGroup[covariateData$methylGroup > 0.75] = 2
    covariateData$methylGroup[(covariateData$methylGroup >= 0.25) & (covariateData$methylGroup <= 0.75)] = NA
    
    kmFit = survfit(patSurv ~ methylGroup, data=covariateData)
    kmPlot = ggsurvplot(kmFit, risk.table = TRUE, conf.int = FALSE)
    # cumcensor = TRUE, cumevents = TRUE,
    # a$plot = a$plot + 
        # scale_color_discrete(name="Strata", labels=c("Low methylation score", 
        #                                                    "High methylation score"), breaks=c("red", "blue")) 
        #theme(legend.text = element_text(c("Low methylation score", "High methylation score")))
    svg(ffPlot(paste0(plotSubdir, "kmPlotTraining", abbrevName[i], ".svg")))
    kmPlot
    dev.off()

}



################################################################################
# part of Panel B
# validation data
vGenomicSignal = methylMat[, -trainDataInd]
vSampleLabels = as.numeric(allSampleLabels[-trainDataInd])
colsToAnnotate = "cancerStage"

valMeta = pMeta[-trainDataInd, ]


# "-" in name causes error
colnames(vGenomicSignal) = gsub(pattern = "-", replacement = "_", x = colnames(vGenomicSignal))

mBySampleDF = aggregateSignalGRList(signal=vGenomicSignal, 
                     signalCoord=signalCoord, GRList=GRList[topRSInd[thisRSInd]], 
                     signalCol = colnames(vGenomicSignal), 
                     scoringMetric = "regionMean", verbose = TRUE)
# convert to long format for plotting
longMByS = transpose(mBySampleDF) 
colnames(longMByS) = abbrevName
longMByS$subjectID = colnames(mBySampleDF)
longMByS = longMByS[1:ncol(vGenomicSignal), ]
longMByS = cbind(longMByS, cancerStage=as.factor(vSampleLabels))
longMByS = gather(data = longMByS, "regionSetName", "methylScore", abbrevName)
mByStagePlot = ggplot(data=longMByS, mapping = aes(x=cancerStage, y=methylScore)) + geom_violin() + facet_wrap("regionSetName") + 
    geom_smooth(mapping = aes(x = as.numeric(cancerStage), y=methylScore)) + xlab("Cancer stage") + ylab("Methylation proportion")
mByStagePlot
ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStageValidation.svg")), plot = mByStagePlot, device = "svg")


# individual RS plots
for (i in seq_along(abbrevName)) {
    mByStagePlotSingle = ggplot(data=filter(longMByS, regionSetName == abbrevName[i]), 
                                mapping = aes(x=cancerStage, y=methylScore)) + 
        geom_violin() +
        geom_smooth(mapping = aes(x = as.numeric(cancerStage), y=methylScore)) + xlab("Cancer stage") + ylab("Methylation proportion")
    mByStagePlotSingle
    ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStageValidation", abbrevName[i], ".svg")), 
           plot = mByStagePlotSingle, device = "svg", width=plotWidth, height = plotHeight / 2,
           units = plotUnit)
}


tmp = rep(-999, nrow(mBySampleDF))
testResDF = data.frame(corPVal=tmp, 
                        corCoef=tmp, 
                        coxPVal=tmp, 
                        coxEffect=tmp)

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
    patSurv = Surv(covariateData$lastDate / (365/12), covariateData$vital_status)
    cModel = coxph(patSurv ~ years_to_birth + gender + meanMethyl + methylScore, data = covariateData)
    print(cModel)
    
    sink(file = ffPlot(paste0(plotSubdir, "coxphModelValidation", abbrevName[i], ".txt")))
    print(cModel)
    print(summary(cModel))
    sink()
    
    # create forest plot for hazard ratios
    fPlot = ggforest(cModel)
    ggsave(filename = ffPlot(paste0(plotSubdir, "forestPlot_validation_", abbrevName[i], ".svg")), 
           plot = fPlot, device = "svg", width = plotWidth *(5/4), height = plotHeight, units = plotUnit)
    
    
    # create groups
    kmThresh = 0.25
    covariateData$methylGroup = ecdf(covariateData$methylScore)(covariateData$methylScore)
    covariateData$methylGroup[covariateData$methylGroup < kmThresh] = 0
    covariateData$methylGroup[covariateData$methylGroup > (1-kmThresh)] = 2
    covariateData$methylGroup[(covariateData$methylGroup >= kmThresh) & (covariateData$methylGroup <= (1-kmThresh))] = NA
    
    kmFit = survfit(patSurv ~ methylGroup, data=covariateData)
    kPlot = ggsurvplot(kmFit, risk.table = TRUE)
    kPlot
    pdf(ffPlot(paste0(plotSubdir, abbrevName[i], "_Kaplan", kmThresh, "Validation", dataID, ".pdf")))
        print(kPlot)
    dev.off()
    # ggsave(filename = ffPlot(paste0(plotSubdir, "EZH2", "Kaplan", kmThresh, dataID, ".svg")) , plot = kPlot, device = "svg")
    
    plot(kmFit)
}


#############################################################################