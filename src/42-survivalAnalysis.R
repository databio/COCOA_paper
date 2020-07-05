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

# # test how much methylation level in each region set is correlated with
# # methylation level in other region sets
# corMat = cor(t(mBySampleDF[, 1:ncol(genomicSignal)]))

# part of Panel B

##########################
# object to save results
myVar = c("methylScore", "age", "sex", 
          "genomeMethyl") # name of variables in my results data.frame
modelVar = c("methylScore", "years_to_birth","gendermale", "meanMethyl") #name of variables in the model 
perModelNum = length(myVar)
rsNum = length(abbrevName)
cancerID = "KIRC"
outerCoxResults = list()
outerSpearmanResults = list()

###########################################################################

# @param hrScaleMod character e.g. "0.01" for methylation levels
# @param ggExpr character ggplot expression beginning with + (e.g. "+ xlab(..)")
jForestPlot = function(plotData, hrScaleMod=NULL, 
                       ggExpr=NULL, minPVal=floor(log10(min(plotData$coxPVal)))) {
    fPlot <- ggplot(data=plotData,
                    aes(x=variableName, y=get(paste0("coxHRMean", hrScaleMod)), 
                        ymin=get(paste0("coxHRLower", hrScaleMod)), ymax=get(paste0("coxHRUpper", hrScaleMod)))) +
        geom_pointrange(aes(col=log10(coxPVal))) + 
        geom_hline(yintercept=1, lty=2) + # facet_wrap(facetVar~., scales = "free") +
        coord_flip() +
        xlab("Variable") + 
        # scale_y_log10() +
        theme_classic() +
        # 
        scale_color_gradient2(low="red", mid="orange", high="gray", midpoint = minPVal/2, limits=c(minPVal, 0))
    # scale_color_manual(values = c("red", "darkorange", "darkgray")) +
    if (!is.null(ggExpr)) {
        fPlot = eval(parse(text = paste0("fPlot", ggExpr)))
    }
    
    
    return(fPlot)
}

#############################################################################
dataSetType = c("training", "testing")
colsToAnnotate = "cancerStage"

a=coxPHPipeline(genomicSignal=methylMat, signalCoord=signalCoord, GRList=GRList[1:2],
                         abbrevName=names(GRList[1:2]), pMeta=pMeta, 
                         covariates=c("years_to_birth", "gender", "methylScore"), 
              covAbbrev= c("age", "gender", "methylScore"), dataID="", 
              plotSubdir=plotSubdir)

# this function 
# @param character. covariates
# @param genomicSignal DNA methylation data
# @param data.frame. pMeta patient metadata. Must have columns: vital_status,
# days_to_last_followup, days_to_death
# @param cancerID e.g. KIRC
# @return list. data.frame with cox variable info, forest plot

coxPHPipeline = function(genomicSignal, signalCoord, GRList,
                         abbrevName=names(GRList), pMeta, 
                         covariates=c("years_to_birth", "gender", "methylScore"), covAbbrev=c(c("age", "gender"), "methylScore"),
                         dataID="", plotSubdir=getwd()) {
    
    
    thisMeta = pMeta
    # methylScore is the variable of interest, created later
    survFormula = as.formula(paste0("patSurv ", "~ ", paste0(c(covariates, "methylScore"), collapse = " + ")))

    # make data.frame to rbind() results to
    coxResults = data.frame(cancerID = NA, 
                            rsName=NA,
                            variableName=NA,
                            obsNum=NA, # number of observations
                            coxCoef=NA,
                            coxCoefSE=NA,
                            coxHRMean=NA, # exp(coxParam)
                            coxHRLower=NA, 
                            coxHRUpper=NA, 
                            coxPVal=NA,
                            schoRho=NA,
                            schoX2=NA,
                            schoPVal=NA,
                            eventNum=NA, # deaths
                            stringsAsFactors = FALSE)
    
    ###########################################################################
    # # making plots of correlation between variable of interest and survival
    # sampleLabels = as.numeric(sampleLabels)
    # spearResults = data.frame(cancerID = rep(cancerID, length(abbrevName)), 
    #                           rsName=rep(NA, length(abbrevName)),
    #                           correlation=rep(NA, length(abbrevName)),
    #                           spearmanPVal=rep(NA, length(abbrevName)),
    #                           stringsAsFactors = FALSE)
    # 
    # # convert to long format for plotting
    # longMByS = transpose(mBySampleDF) 
    # colnames(longMByS) = abbrevName
    # longMByS$subjectID = colnames(mBySampleDF)
    # longMByS = longMByS[1:ncol(genomicSignal), ]
    # longMByS = cbind(longMByS, cancerStage=as.factor(sampleLabels))
    # longMByS = gather(data = longMByS, "regionSetName", "methylScore", abbrevName)
    # mByStagePlot = ggplot(data=longMByS, mapping = aes(x=cancerStage, y=methylScore)) + geom_violin() + facet_wrap("regionSetName") + 
    #     geom_smooth(mapping = aes(x = as.numeric(cancerStage), y=methylScore)) + xlab("Cancer stage") + ylab("Methylation proportion")
    # mByStagePlot
    # ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStage", ".svg")), plot = mByStagePlot, device = "svg")
    # 
    # # individual RS plots
    # for (i in seq_along(abbrevName)) {
    #     mByStagePlotSingle = ggplot(data=filter(longMByS, regionSetName == abbrevName[i]), 
    #                                 mapping = aes(x=cancerStage, y=methylScore)) + 
    #         geom_violin() +
    #         geom_smooth(mapping = aes(x = as.numeric(cancerStage), y=methylScore)) + xlab("Cancer stage") + ylab("Methylation proportion")
    #     mByStagePlotSingle
    #     ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStage", "_", abbrevName[i], ".svg")), 
    #            plot = mByStagePlotSingle, device = "svg", width=plotWidth, height = plotHeight / 2,
    #            units = plotUnit)
    # }
    # 
    # plot(sampleLabels, mBySample)
    # thisSpear = cor.test(x = sampleLabels, mBySample, method="spearman")
    # spearResults[i, "rsName"] = abbrevName[i]
    # spearResults[i, "correlation"] = thisSpear$estimate
    # spearResults[i, "spearmanPVal"] =  thisSpear$p.value
    # 
    # print(cor.test(x = thisMeta$years_to_birth, mBySample))
    # 
    # aliveMethyl = mBySample[thisMeta$vital_status == 0]
    # deadMethyl = mBySample[thisMeta$vital_status == 1]
    # print(wilcox.test(aliveMethyl, deadMethyl, conf.int = TRUE))
    # thisMeta$meanMethyl = colMeans(genomicSignal)
    # 
    
    ############################################################################
    
    ###################################
    # getting average methylation in each region set
    # "-" in name causes error
    colnames(genomicSignal) = gsub(pattern = "-", replacement = "_", x = colnames(genomicSignal))
    
    mBySampleDF = aggregateSignalGRList(signal=genomicSignal, 
                                        signalCoord=signalCoord, GRList=GRList, 
                                        signalCol = colnames(genomicSignal), 
                                        scoringMetric = "regionMean", verbose = TRUE)
    
    # once for each region set
    for (i in 1:nrow(mBySampleDF)) {
        mBySample = as.numeric(mBySampleDF[i, 1:ncol(genomicSignal)])
        thisMeta$methylScore = mBySample
        
        # test whether DNA methylation is associated with survival
        ###### cox proportional hazards model 
        thisMeta$lastDate = thisMeta$days_to_last_followup
        thisMeta$lastDate[is.na(thisMeta$lastDate)] = thisMeta$days_to_death[is.na(thisMeta$lastDate)]
        
        # covariates
        covariateData = thisMeta
        patSurv = Surv(covariateData$lastDate / (365/12), covariateData$vital_status)
        covariateData$patSurv = patSurv
        cModel = coxph(survFormula, data = covariateData)
        print(cModel)
        
        sink(file = ffPlot(paste0(plotSubdir, "coxphModel_", dataID, "_", abbrevName[i], ".txt")))
        print(cModel)
        print(summary(cModel))
        sink()
        
        # create forest plot for hazard ratios
        fPlot = ggforest(cModel, data = covariateData)
        ggsave(filename = ffPlot(paste0(plotSubdir, "forestPlot_", dataID, "_", abbrevName[i], ".svg")), 
               plot = fPlot, device = "svg", width = plotWidth *(5/4), height = plotHeight, units = plotUnit)
        
        modelVar = names(cModel$coefficients)
        myVar = covAbbrev
        thisCoxDF = makeCoxDF(cModel, 
                              modelVar=modelVar, 
                              dfVar=myVar, 
                              returnGlobalStats = TRUE)
        
        thisCoxDF = cbind(cancerID=rep(cancerID, nrow(thisCoxDF)),
                          rsName=rep(abbrevName[i], nrow(thisCoxDF)), 
                          thisCoxDF)
        coxResults = rbind(coxResults, thisCoxDF)
        
        ############################################################
        # # kaplan meier
        # # create groups
        # kmThresh = 0.25
        # covariateData$methylGroup = ecdf(covariateData$methylScore)(covariateData$methylScore)
        # covariateData$methylGroup[covariateData$methylGroup < kmThresh] = 0
        # covariateData$methylGroup[covariateData$methylGroup > (1-kmThresh)] = 2
        # covariateData$methylGroup[(covariateData$methylGroup >= kmThresh) & (covariateData$methylGroup <= (1-kmThresh))] = NA
        # 
        # kmFit = survfit(patSurv ~ methylGroup, data=covariateData)
        # kPlot = ggsurvplot(kmFit, data = covariateData, risk.table = TRUE, pval = TRUE, pval.method = TRUE)
        # kPlot
        # pdf(ffPlot(paste0(plotSubdir, abbrevName[i], "_Kaplan", kmThresh, "_", dataID, ".pdf")))
        # print(kPlot)
        # dev.off()
        # # ggsave(filename = ffPlot(paste0(plotSubdir, "EZH2", "Kaplan", kmThresh, dataID, ".svg")) , plot = kPlot, device = "svg")
        # 
        # plot(kmFit)
        ##############################################################
    }
    coxResults = coxResults[-1, ]
    
    # convert hazard ratios from 1 scale to 0.01 scale
    coxResults$coxHRMean0.01=coxResults$coxHRMean ^ (1/100)
    coxResults$coxHRUpper0.01=coxResults$coxHRUpper ^ (1/100)
    coxResults$coxHRLower0.01=coxResults$coxHRLower ^ (1/100)
    
    # correct by variable
    coxResults = data.table(coxResults)
    coxResults[, holmCoxPVal := p.adjust(coxPVal, method = "BH"), by=variableName]
    coxResults = as.data.frame(coxResults)
    all(filter(coxResults, variableName == "sex")$holmCoxPVal == p.adjust(filter(coxResults, variableName == "sex")$coxPVal, "BH"))
    #summarise(a, test=p.adjust(coxPVal, method = "BH"))
    
    coxResults$sigType = rep(x = "p > 0.05", nrow(coxResults))
    # coxResults$sigType[coxResults$coxPVal < 0.05] = "Uncorrected p < 0.05"
    coxResults$sigType[coxResults$holmCoxPVal < 0.05] = "Corrected p < 0.05"
    coxResults$sigType = factor(coxResults$sigType, 
                                levels = c("Corrected p < 0.05", 
                                           "Uncorrected p < 0.05", 
                                           "p > 0.05"))
    
    ###########################################################################
    
    for (i in seq_along(abbrevName)) {
        
        plotData = filter(coxResults, (rsName == abbrevName[i])) 
        
        # create dummy variable to use for faceting and making different scales
        plotData$facetVar = rep("normal scale", nrow(plotData))
        plotData$facetVar[plotData$variableName %in% c("genomeMethyl", "methylScore")] = 
            "special scale"
        # for color scale
        minPVal = floor(log10(min(plotData$coxPVal)))
        
        fPlot = jForestPlot(filter(plotData, variableName %in% c("genomeMethyl", "methylScore")), 
                            hrScaleMod="0.01", minPVal = minPVal) +
            #ggExpr = paste0("+ scale_color_gradient2(low='red', mid='orange', high='gray', midpoint = min(log10(plotData$coxPVal))/2, limits=c(0,", minPVal,"))")) + 
            ylab("Hazard ratio per 0.01 change in DNA methylation level (95% CI)") +
            scale_y_continuous(breaks=seq(from=0.8, to=2, by=0.2)) + theme(aspect.ratio = 1)
        
        if (nrow(filter(plotData, (as.character(sigType) == "Corrected p < 0.05") &
                        (variableName %in% c("genomeMethyl", "methylScore")))) > 0) {
            fPlot = fPlot + geom_text(mapping = aes(x=variableName, y= 1.8),
                                      data = filter(plotData, (as.character(sigType) == "Corrected p < 0.05") &
                                                        (variableName %in% c("genomeMethyl", "methylScore"))), 
                                      label="*", size = 6)
        }
        
        
        fPlot
        ggsave(filename = ffPlot(paste0(plotSubdir, "coxHazardRatios1Val_", dataID, "_", cancerID, "_", abbrevName[i], ".svg")), 
               plot = fPlot, device = "svg", height = 50, width=100, units = "mm")
        
        fPlot = jForestPlot(filter(plotData, variableName %in% c("sex", "age")), minPVal=minPVal) +
            # ggExpr = paste0("+ scale_color_gradient2(low='red', mid='orange', high='gray', midpoint = min(log10(plotData$coxPVal))/2, limits=c(0,", minPVal,"))")) + 
            ylab("Hazard ratio") + theme(aspect.ratio = 1)
        
        if (nrow(filter(plotData, (variableName %in% c("sex", "age")) &
                        (as.character(sigType) == "Corrected p < 0.05"))) > 0) {
            fPlot = fPlot + geom_text(mapping = aes(x=variableName, y= 1.2),
                                      data = filter(plotData, (variableName %in% c("sex", "age")) &
                                                        (as.character(sigType) == "Corrected p < 0.05")),
                                      label="*", size = 6) #+
            
        }
        #scale_y_continuous(breaks=seq(from=0.8, to=1.1, by=0.2))
        fPlot
        ggsave(filename = ffPlot(paste0(plotSubdir, "coxHazardRatios2Val_", dataID, "_", cancerID, "_", 
                                        abbrevName[i], ".svg")), 
               plot = fPlot, device = "svg", height = 50, width=100, units = "mm")
    }
    
    return(list(coxResults))
}

###################################################################################

for (dataInd in seq_along(dataSetType)) {
    
    # make data.frame to rbind() results to
    coxResults = data.frame(dataSetType=NA, 
                            cancerID = NA, 
                            rsName=NA,
                            variableName=NA,
                            obsNum=NA, # number of observations
                            coxCoef=NA,
                            coxCoefSE=NA,
                            coxHRMean=NA, # exp(coxParam)
                            coxHRLower=NA, 
                            coxHRUpper=NA, 
                            coxPVal=NA,
                            schoRho=NA,
                            schoX2=NA,
                            schoPVal=NA,
                            eventNum=NA, # deaths
                            stringsAsFactors = FALSE)
    spearResults = data.frame(dataSetType=rep(dataSetType[dataInd], length(abbrevName)), 
                              cancerID = rep(cancerID, length(abbrevName)), 
                              rsName=rep(NA, length(abbrevName)),
                              correlation=rep(NA, length(abbrevName)),
                              spearmanPVal=rep(NA, length(abbrevName)),
                              stringsAsFactors = FALSE)
    
    ###################################
    # getting average methylation in each region set
    
    if (dataSetType[dataInd] == "training") {
        # training data
        genomicSignal = methylMat[, trainDataInd]
        sampleLabels = as.numeric(allSampleLabels[trainDataInd])
        thisMeta = pMeta[trainDataInd, ]
    } else {
        # validation data
        genomicSignal = methylMat[, -trainDataInd]
        sampleLabels = as.numeric(allSampleLabels[-trainDataInd])
        thisMeta = pMeta[-trainDataInd, ]
    }

    
    ###################################
    # "-" in name causes error
    colnames(genomicSignal) = gsub(pattern = "-", replacement = "_", x = colnames(genomicSignal))
    
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
    ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStage", dataSetType[dataInd], "_", ".svg")), plot = mByStagePlot, device = "svg")
    
    
    # individual RS plots
    for (i in seq_along(abbrevName)) {
        mByStagePlotSingle = ggplot(data=filter(longMByS, regionSetName == abbrevName[i]), 
                                    mapping = aes(x=cancerStage, y=methylScore)) + 
            geom_violin() +
            geom_smooth(mapping = aes(x = as.numeric(cancerStage), y=methylScore)) + xlab("Cancer stage") + ylab("Methylation proportion")
        mByStagePlotSingle
        ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStage", dataSetType[dataInd], "_", abbrevName[i], ".svg")), 
               plot = mByStagePlotSingle, device = "svg", width=plotWidth, height = plotHeight / 2,
               units = plotUnit)
    }
    
    #############################################
    # once for each region set
    for (i in 1:nrow(mBySampleDF)) {
        mBySample = as.numeric(mBySampleDF[i, 1:ncol(genomicSignal)])
        thisMeta$methylScore = mBySample
        plot(sampleLabels, mBySample)
        thisSpear = cor.test(x = sampleLabels, mBySample, method="spearman")
        spearResults[i, "rsName"] = abbrevName[i]
        spearResults[i, "correlation"] = thisSpear$estimate
        spearResults[i, "spearmanPVal"] =  thisSpear$p.value
        outerSpearmanResults[[dataInd]] = spearResults
        
        print(cor.test(x = thisMeta$years_to_birth, mBySample))
        
        aliveMethyl = mBySample[thisMeta$vital_status == 0]
        deadMethyl = mBySample[thisMeta$vital_status == 1]
        print(wilcox.test(aliveMethyl, deadMethyl, conf.int = TRUE))
        thisMeta$meanMethyl = colMeans(genomicSignal)
        
        
        # test whether DNA methylation is associated with survival
        ###### cox proportional hazards model 
        thisMeta$lastDate = thisMeta$days_to_last_followup
        thisMeta$lastDate[is.na(thisMeta$lastDate)] = thisMeta$days_to_death[is.na(thisMeta$lastDate)]
        
        # covariates
        covariateData = thisMeta
        patSurv = Surv(covariateData$lastDate / (365/12), covariateData$vital_status)
        cModel = coxph(patSurv ~ years_to_birth + gender + meanMethyl + methylScore, data = covariateData)
        print(cModel)
        
        sink(file = ffPlot(paste0(plotSubdir, "coxphModel_", dataSetType[dataInd], "_", abbrevName[i], ".txt")))
        print(cModel)
        print(summary(cModel))
        sink()
        
        # create forest plot for hazard ratios
        fPlot = ggforest(cModel)
        ggsave(filename = ffPlot(paste0(plotSubdir, "forestPlot_", dataSetType[dataInd], "_", abbrevName[i], ".svg")), 
               plot = fPlot, device = "svg", width = plotWidth *(5/4), height = plotHeight, units = plotUnit)
        
        modelVar = c("years_to_birth", "gendermale", "meanMethyl", "methylScore")
        myVar = c("age", "sex", "genomeMethyl", "methylScore")
        thisCoxDF = makeCoxDF(cModel, 
                              modelVar=modelVar, 
                              dfVar=myVar, 
                              returnGlobalStats = TRUE)
        
        thisCoxDF = cbind(dataSetType=rep(dataSetType[dataInd], nrow(thisCoxDF)), 
                          cancerID=rep(cancerID, nrow(thisCoxDF)),
                          rsName=rep(abbrevName[i], nrow(thisCoxDF)), 
                          thisCoxDF)
        coxResults = rbind(coxResults, thisCoxDF)
        
        
        # create groups
        kmThresh = 0.25
        covariateData$methylGroup = ecdf(covariateData$methylScore)(covariateData$methylScore)
        covariateData$methylGroup[covariateData$methylGroup < kmThresh] = 0
        covariateData$methylGroup[covariateData$methylGroup > (1-kmThresh)] = 2
        covariateData$methylGroup[(covariateData$methylGroup >= kmThresh) & (covariateData$methylGroup <= (1-kmThresh))] = NA
        
        kmFit = survfit(patSurv ~ methylGroup, data=covariateData)
        kPlot = ggsurvplot(kmFit, risk.table = TRUE, pval = TRUE, pval.method = TRUE)
        kPlot
        pdf(ffPlot(paste0(plotSubdir, abbrevName[i], "_Kaplan", kmThresh, dataSetType[dataInd], "_", dataID, ".pdf")))
        print(kPlot)
        dev.off()
        # ggsave(filename = ffPlot(paste0(plotSubdir, "EZH2", "Kaplan", kmThresh, dataID, ".svg")) , plot = kPlot, device = "svg")
        
        plot(kmFit)
    }
    
        coxResults = coxResults[-1, ]
    
    # convert hazard ratios from 1 scale to 0.01 scale
    coxResults$coxHRMean0.01=coxResults$coxHRMean ^ (1/100)
    coxResults$coxHRUpper0.01=coxResults$coxHRUpper ^ (1/100)
    coxResults$coxHRLower0.01=coxResults$coxHRLower ^ (1/100)
    
    # correct by variable
    
    coxResults = data.table(coxResults)
    coxResults[, holmCoxPVal := p.adjust(coxPVal, method = "BH"), by=variableName]
    coxResults = as.data.frame(coxResults)
    all(filter(coxResults, variableName == "sex")$holmCoxPVal == p.adjust(filter(coxResults, variableName == "sex")$coxPVal, "BH"))
    #summarise(a, test=p.adjust(coxPVal, method = "BH"))
    
    coxResults$sigType = rep(x = "p > 0.05", nrow(coxResults))
    # coxResults$sigType[coxResults$coxPVal < 0.05] = "Uncorrected p < 0.05"
    coxResults$sigType[coxResults$holmCoxPVal < 0.05] = "Corrected p < 0.05"
    coxResults$sigType = factor(coxResults$sigType, 
                                levels = c("Corrected p < 0.05", 
                                           "Uncorrected p < 0.05", 
                                           "p > 0.05"))
    
    outerCoxResults[[dataInd]] = coxResults
    
    ###########################################################################
    # make forest plots

    for (i in seq_along(abbrevName)) {
        
        
        plotData = filter(coxResults, (rsName == abbrevName[i])) 
        
        # create dummy variable to use for faceting and making different scales
        plotData$facetVar = rep("normal scale", nrow(plotData))
        plotData$facetVar[plotData$variableName %in% c("genomeMethyl", "methylScore")] = 
            "special scale"
        # for color scale
        minPVal = floor(log10(min(plotData$coxPVal)))
        
        fPlot = jForestPlot(filter(plotData, variableName %in% c("genomeMethyl", "methylScore")), 
                            hrScaleMod="0.01", minPVal = minPVal) +
                            #ggExpr = paste0("+ scale_color_gradient2(low='red', mid='orange', high='gray', midpoint = min(log10(plotData$coxPVal))/2, limits=c(0,", minPVal,"))")) + 
            ylab("Hazard ratio per 0.01 change in DNA methylation level (95% CI)") +
            scale_y_continuous(breaks=seq(from=0.8, to=2, by=0.2)) + theme(aspect.ratio = 1)
            
        if (nrow(filter(plotData, (as.character(sigType) == "Corrected p < 0.05") &
                        (variableName %in% c("genomeMethyl", "methylScore")))) > 0) {
            fPlot = fPlot + geom_text(mapping = aes(x=variableName, y= 1.8),
                                     data = filter(plotData, (as.character(sigType) == "Corrected p < 0.05") &
                                                       (variableName %in% c("genomeMethyl", "methylScore"))), 
                                     label="*", size = 6)
        }
        
        
        fPlot
        ggsave(filename = ffPlot(paste0(plotSubdir, "coxHazardRatios1Val_", dataSetType[dataInd], "_", cancerID, "_", abbrevName[i], ".svg")), 
               plot = fPlot, device = "svg", height = 50, width=100, units = "mm")
        
        fPlot = jForestPlot(filter(plotData, variableName %in% c("sex", "age")), minPVal=minPVal) +
                            # ggExpr = paste0("+ scale_color_gradient2(low='red', mid='orange', high='gray', midpoint = min(log10(plotData$coxPVal))/2, limits=c(0,", minPVal,"))")) + 
        ylab("Hazard ratio") + theme(aspect.ratio = 1)
        
        if (nrow(filter(plotData, (variableName %in% c("sex", "age")) &
                        (as.character(sigType) == "Corrected p < 0.05"))) > 0) {
            fPlot = fPlot + geom_text(mapping = aes(x=variableName, y= 1.2),
                                      data = filter(plotData, (variableName %in% c("sex", "age")) &
                                                        (as.character(sigType) == "Corrected p < 0.05")),
                                      label="*", size = 6) #+
            
        }
        #scale_y_continuous(breaks=seq(from=0.8, to=1.1, by=0.2))
        fPlot
        ggsave(filename = ffPlot(paste0(plotSubdir, "coxHazardRatios2Val_", dataSetType[dataInd], "_", cancerID, "_", 
                                        abbrevName[i], ".svg")), 
               plot = fPlot, device = "svg", height = 50, width=100, units = "mm")
    }
    
    #############################################################################
}

coxResults = rbindlist(outerCoxResults)
write.csv(x = coxResults, file = ffSheets("KIRC_coxResults.csv"), quote = FALSE)

spearmanResults = rbindlist(outerSpearmanResults)
write.csv(x = spearmanResults, file = ffSheets("KIRC_spearmanCancerStageResults.csv"), quote = FALSE)
