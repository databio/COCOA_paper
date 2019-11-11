# look at methylation in EZH2 binding regions in many different cancers

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(curatedTCGAData)
library(TCGAutils)

library(ExperimentHub)
ExperimentHub::getExperimentHubOption("CACHE")
setExperimentHubOption("CACHE", "/scratch/jtl2hk/.cache/ExperimentHub")

plotWidth = 100
plotHeight = 100
plotUnit = "mm"

#############################################################################
loadGRList(genomeV = "hg19")

myTopRSNames = c("wgEncodeAwgTfbsBroadHsmmtEzh239875UniPk.narrowPeak",
"wgEncodeAwgTfbsSydhNt2d1Suz12UcdUniPk.narrowPeak")
myTopRS = GRList[myTopRSNames]
abbrevNames = c("EZH2", "SUZ12")
##############################################################################

# do for each cancer type
cancerID= c("ACC", "BLCA", "BRCA", "CESC", "CHOL", 
            "COAD", "DLBC", "ESCA", "GBM", "HNSC", 
            "KICH", "KIRC", "KIRP", "LAML", "LGG", 
            "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG",
"PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") 

spearCor=rep(NA, length(cancerID)*length(myTopRS))
spearCorDF = data.frame(cancerID = rep(cancerID, each=length(myTopRS)), 
                        rsName=spearCor, spearCor, pVal=spearCor, 
                        stringsAsFactors = FALSE)
for (i in seq_along(cancerID)) {
    
    allSampleLabels = NULL
    # assigns methylMat, signalCoord, pMeta, allSampleLabels to environment
    fxCode = loadProcessTCGAMethyl(cancerID[i])
    
    if (is.null(fxCode)) {
        next()
    }
    if (!all(c("days_to_last_followup", "days_to_death", 
               "years_to_birth", "gender", "vital_status") %in% colnames(pMeta))) {
        next()
    }

    
    genomicSignal = methylMat
    
    # get average methylation in EZH2 and SUZ12 regions
    # "-" in name causes error
    colnames(genomicSignal) = gsub(pattern = "-", replacement = "_", x = colnames(genomicSignal))
    # average methylation per sample in these regions
    mBySampleDF = aggregateSignalGRList(signal=genomicSignal, 
                                        signalCoord=signalCoord, GRList=myTopRS, 
                                        signalCol = colnames(genomicSignal), 
                                        scoringMetric = "regionMean", verbose = TRUE)
    # convert to long format for plotting
    longMByS = transpose(mBySampleDF) 
    colnames(longMByS) = abbrevName
    longMByS$subjectID = colnames(mBySampleDF)
    longMByS = longMByS[1:ncol(genomicSignal), ]
    if (!is.null(allSampleLabels)) {
        longMByS = cbind(longMByS, cancerStage=as.factor(allSampleLabels))
    }
    
    longMByS = gather(data = longMByS, "regionSetName", "methylScore", abbrevName)
    if (!is.null(allSampleLabels)) {
        mByStagePlot = ggplot(data=longMByS, mapping = aes(x=cancerStage, y=methylScore)) + geom_violin() + facet_wrap("regionSetName") + 
            geom_smooth(mapping = aes(x = as.numeric(cancerStage), y=methylScore)) + xlab("Cancer stage") + ylab("Methylation proportion")
        mByStagePlot
        ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStage_", cancerID[i], ".svg")), 
               plot = mByStagePlot, device = "svg")
        
        # individual RS plots and spearman correlation
        for (j in seq_along(abbrevName)) {
            mByStagePlotSingle = ggplot(data=filter(longMByS, regionSetName == abbrevName[j]), 
                                        mapping = aes(x=cancerStage, y=methylScore)) + 
                geom_violin() +
                geom_smooth(mapping = aes(x = as.numeric(cancerStage), y=methylScore)) + xlab("Cancer stage") + ylab("Methylation proportion")
            mByStagePlotSingle
            ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStage", abbrevName[j],"_", cancerID[i], ".svg")), 
                   plot = mByStagePlotSingle, device = "svg", width=plotWidth, height = plotHeight / 2,
                   units = plotUnit)
            
            corRes = cor.test(x = allSampleLabels, 
                              as.vector(mBySample[j, 1:ncol(genomicSignal)]), method = "spearman")
            spearCorDF[2*(i-1)+j, "rsName"] = abbrevName[j]
            spearCorDF[2*(i-1)+j, "spearCor"] = corRes$estimate
            spearCorDF[2*(i-1)+j, "pVal"] = corRes$p.value
        }
    }
    
    
    # test association between avg methyl. and cancer stage, and survival
    
    # data.frame to store results 
    tmp = rep(-999, nrow(mBySampleDF))
    trainResDF = data.frame(corPVal=tmp, 
                            corCoef=tmp, 
                            coxPVal=tmp, 
                            coxEffect=tmp, 
                            coxModel=tmp)
    
    # once for each region set
    for (j in 1:nrow(mBySampleDF)) {
        mBySample = as.numeric(mBySampleDF[j, 1:ncol(genomicSignal)])
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
        
        sink(file = ffPlot(paste0(plotSubdir, "coxphModel", abbrevName[j], "_", cancerID[i], ".txt")))
        print(myModel)
        print(summary(myModel))
        sink()
        
        
        # kaplan meier
        # create groups
        covariateData$methylGroup = ecdf(covariateData$methylScore)(covariateData$methylScore)
        covariateData$methylGroup[covariateData$methylGroup < 0.3333] = 0
        covariateData$methylGroup[covariateData$methylGroup > 0.6666] = 2
        covariateData$methylGroup[(covariateData$methylGroup >= 0.3333) & (covariateData$methylGroup <= 0.6666)] = NA
        
        kmFit = survfit(patSurv ~ methylGroup, data=covariateData)
        kmPlot = ggsurvplot(kmFit, risk.table = TRUE, conf.int = TRUE)
        # cumcensor = TRUE, cumevents = TRUE,
        # a$plot = a$plot + 
        # scale_color_discrete(name="Strata", labels=c("Low methylation score", 
        #                                                    "High methylation score"), breaks=c("red", "blue")) 
        #theme(legend.text = element_text(c("Low methylation score", "High methylation score")))
        svg(ffPlot(paste0(plotSubdir, "kmPlot", abbrevName[j], "_", cancerID[i], ".svg")))
        kmPlot
        dev.off()
        
    }
    
    # store in DF

    
}
write.csv(spearCorDF, file = ffSheets("EZH2_TCGA_cancer_stage.csv"), quote = FALSE, 
          row.names = FALSE, col.names = TRUE)

