# look at methylation in EZH2 binding regions in many different cancers

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(curatedTCGAData)
library(TCGAutils)
library(tidyr)
library(survival)
library(survminer)

library(ExperimentHub)
ExperimentHub::getExperimentHubOption("CACHE")
setExperimentHubOption("CACHE", "/scratch/jtl2hk/.cache/ExperimentHub")

plotWidth = 100
plotHeight = 100
plotUnit = "mm"

plotSubdir = "99-ezh2Analysis/"
if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

#############################################################################
loadGRList(genomeV = "hg19")

##############################################################################
# create consensus polycomb region set 

# combine top region sets from KIRC analysis
dataID = "kircMethyl214"
variationMetric = "spearmanCor" 
simpleCache(paste0("rsScores_", dataID, "_", variationMetric), assignToVariable = "rsScores")

# ezh2 and suz12
topEzh2 = arrange(rsScores, desc(cancerStage))$rsName[1:11]
topEzh2 = topEzh2[grepl(pattern = "ezh2|suz12", x = topEzh2, ignore.case = TRUE)]
topEzh2GRList = GRList[topEzh2]
mergedTopEzh2 = unlist(topEzh2GRList, recursive = TRUE, use.names = TRUE)
mergedTopEzh2 = reduce(mergedTopEzh2)

# myTopRSNames = c("wgEncodeAwgTfbsBroadHsmmtEzh239875UniPk.narrowPeak",
#                  "wgEncodeAwgTfbsSydhNt2d1Suz12UcdUniPk.narrowPeak")
myTopRS = mergedTopEzh2
# myTopRS must be a GRList for following code
myTopRS = GRangesList(myTopRS)
abbrevName = c("polycomb")
length(mergedTopEzh2) # 10994

###############################################################################

# @param dfVar character. The column names of the new df

makeCoxDF <- function(cModel, modelVar, dfVar, returnGlobalStats=TRUE) {
    
    if (returnGlobalStats) {
        modelVar = c(modelVar, "GLOBAL")
        dfVar = c(dfVar, "globalStats")
    }
    perModelNum = length(modelVar)
    
    phTest = cox.zph(cModel)
    
    coxVar=rep(NA, length(dfVar))
    # cancerID = rep(cancerID, each=perModelNum), rsName=coxVar,
    coxDF = data.frame(variableName=dfVar,
                       obsNum=coxVar, # number of observations
                       coxHRMean=coxVar, 
                       coxHRLower=coxVar, 
                       coxHRUpper=coxVar, 
                       coxPVal=coxVar,
                       schoRho=coxVar,
                       schoX2=coxVar,
                       schoPVal=coxVar,
                       eventNum=coxVar, # deaths
                       stringsAsFactors = FALSE)
    
    for (i in seq_along(modelVar)) {
        if (modelVar[i] != "GLOBAL") {
            coxDF[i, "variableName"] = dfVar[i]
            # coxDF[i, "obsNum"] =
            coxDF[i, "coxHRMean"] = summary(cModel)$coefficients[modelVar[i], "exp(coef)"]
            coxDF[i, "coxHRUpper"] = summary(cModel)$conf.int[modelVar[i], "upper .95"]
            coxDF[i, "coxHRLower"] = summary(cModel)$conf.int[modelVar[i], "lower .95"]
            coxDF[i, "coxPVal"] = summary(cModel)$coefficients[modelVar[i], "Pr(>|z|)"]
            # coxDF[i, "schoRho"] = phTest$table[modelVar[i],"rho"]
            coxDF[i, "schoX2"] = phTest$table[modelVar[i],"chisq"]
            coxDF[i, "schoPVal"] = phTest$table[modelVar[i],"p"]
        } else {
            coxDF[i, "variableName"] = dfVar[i]
            # log ratio/likelihood ratio test
            coxDF[i, "coxPVal"] = summary(cModel)$logtest["pvalue"]  
            coxDF[i, "schoX2"] = phTest$table[modelVar[i],"chisq"]
            coxDF[i, "schoPVal"] = phTest$table[modelVar[i],"p"]
            coxDF[i, "eventNum"] = summary(cModel)$nevent # deaths
        }
        
        
    }
    
    return(coxDF)
}

###############################################################################

# do for each cancer type
cancerID= c("ACC", "BLCA", "BRCA", "CESC", "CHOL", 
            "COAD", "DLBC", "ESCA", "GBM", "HNSC", 
            "KICH", "KIRC", "KIRP", "LAML", "LGG", 
            "LIHC", "LUAD", "LUSC", "MESO", "OV", 
            "PAAD", "PCPG", "PRAD", "READ", "SARC", 
            "SKCM", "STAD", "TGCT", "THCA", "THYM", 
            "UCEC", "UCS", "UVM") 

spearCor=rep(NA, length(cancerID)*length(myTopRS))
spearCorDF = data.frame(cancerID = rep(cancerID, each=length(myTopRS)), 
                        rsName=spearCor, spearCor, spearmanPVal=spearCor, 
                        coxExp=spearCor, coxPVal=spearCor,
                        stringsAsFactors = FALSE)

# number of cancers * number of region sets * (cox variables + 1 for whole model stats)
# number of cox variables is 4 (sex, genome-wide average methylation, age, methylScore) + 1 = 5
# global is not a variable, just a place to record stats about the model as a whole
# if 
# redefined inside loop
myVar = c("methylScore", "age", "sex", 
          "genomeMethyl") # name of variables in my results data.frame
modelVar = c("methylScore", "years_to_birth","gender", "meanMethyl") #name of variables in the model 
perModelNum = length(myVar)
rsNum = length(myTopRS)
coxVar=rep(NA, length(cancerID)*length(myTopRS) * perModelNum)
# make data.frame to rbind() results to
coxResults = data.frame(cancerID = NA, 
                        rsName=NA,
                        variableName=NA,
                        obsNum=NA, # number of observations
                        coxHRMean=NA, 
                        coxHRLower=NA, 
                        coxHRUpper=NA, 
                        coxPVal=NA,
                        schoRho=NA,
                        schoX2=NA,
                        schoPVal=NA,
                        eventNum=NA, # deaths
                        stringsAsFactors = FALSE)

# get overlap between myTopRS (based on CpGs covered in epigenetic data)
getOverlap=FALSE
for (i in seq_along(cancerID)) {
    message(i)
    myVar = c("methylScore", "age", "sex", 
              "genomeMethyl") # name of variables in my results data.frame
    modelVar = c("methylScore", "years_to_birth","gender", "meanMethyl") 
    
    allSampleLabels = NULL
    # assigns methylMat, signalCoord, pMeta, allSampleLabels to environment
    fxCode = loadProcessTCGAMethyl(cancerID[i])
    
    if (is.null(fxCode)) {
        next()
    }
    if (!all(c("days_to_last_followup", "days_to_death", 
              "vital_status") %in% colnames(pMeta))) {
        next()
    }
    
    genomicSignal = methylMat
    # higher number should be dead, lower number should be alive (either 0/1 or 1/2)
    if (is(pMeta$vital_status, "character")) {
        pMeta$vital_status = factor(pMeta$vital_status, levels = c("alive", "dead"))
    }

    # if (!is(pMeta$vital_status, "integer")) {
    #     message(paste0(class(pMeta$vital_status), "_", unique(pMeta$vital_status)))
    # }
    
    if (is(pMeta$gender, "character")) {
        pMeta$gender = as.numeric(factor(pMeta$gender, levels = c("female", "male")))
    }
    if (any(is.na(pMeta$vital_status))) {
        naInd = which(is.na(pMeta$vital_status))
        pMeta = pMeta[-naInd, ]
        genomicSignal = genomicSignal[, -naInd]
        if (!is.null(allSampleLabels)) {
            allSampleLabels = allSampleLabels[-naInd]    
        }
        
    }
    
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
        ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStage_", cancerID[i], ".pdf")), 
               plot = mByStagePlot, device = "pdf")
        
        # individual RS plots and spearman correlation
        for (j in seq_along(abbrevName)) {
            mByStagePlotSingle = ggplot(data=filter(longMByS, regionSetName == abbrevName[j]), 
                                        mapping = aes(x=cancerStage, y=methylScore)) + 
                geom_violin() +
                geom_smooth(mapping = aes(x = as.numeric(cancerStage), y=methylScore)) + xlab("Cancer stage") + ylab("Methylation proportion")
            mByStagePlotSingle
            ggsave(filename = ffPlot(paste0(plotSubdir, "methylByStage", abbrevName[j],"_", cancerID[i], ".pdf")), 
                   plot = mByStagePlotSingle, device = "pdf", width=plotWidth, height = plotHeight / 2,
                   units = plotUnit)
            
            corRes = cor.test(x = as.numeric(allSampleLabels), 
                              as.numeric(mBySampleDF[j, 1:ncol(genomicSignal)]), method = "spearman")
            spearCorDF[length(myTopRS)*(i-1)+j, "rsName"] = abbrevName[j]
            spearCorDF[length(myTopRS)*(i-1)+j, "spearCor"] = corRes$estimate
            spearCorDF[length(myTopRS)*(i-1)+j, "spearmanPVal"] = corRes$p.value
        }
    }
    
    
    # test association between avg methyl. and cancer stage, and survival
    
    # # data.frame to store results 
    # spearCorDF


    # these patients have no survival data and so are not helpful for survival
    # analysis despite still being helpful for cancer stage analysis earlier
    removeInd = pMeta$days_to_last_followup == 0
    removeInd[is.na(removeInd)] = FALSE # keep NA rows
    pMeta = pMeta[!removeInd, ]
    if (sum(removeInd) > 0) {
        mBySampleDF = mBySampleDF[, -which(removeInd)]
    }
    genomicSignal = genomicSignal[, !removeInd]
    
    if ("dead" %in% pMeta$vital_status) {
        # dead is 1, alive is 0
        pMeta$vital_status = pMeta$vital_status == "dead"
    }

    
    # once for each region set
    for (j in 1:nrow(mBySampleDF)) {
        mBySample = as.numeric(mBySampleDF[j, 1:ncol(genomicSignal)])
        pMeta$methylScore = mBySample
        pMeta$meanMethyl = colMeans(genomicSignal)

        # aliveMethyl = mBySample[pMeta$vital_status == 0]
        # deadMethyl = mBySample[pMeta$vital_status == 1]
        # print(wilcox.test(aliveMethyl, deadMethyl, conf.int = TRUE))

        # test whether DNA methylation is associated with survival
        ###### cox proportional hazards model
        pMeta$lastDate = pMeta$days_to_last_followup
        pMeta$lastDate[is.na(pMeta$lastDate)] = pMeta$days_to_death[is.na(pMeta$lastDate)]

        # covariates
        covariateData = pMeta
        patSurv = Surv(covariateData$lastDate / (365/12), event=covariateData$vital_status)
        try({
            
            # create model with only variables that are present
            myVar = myVar[modelVar %in% colnames(covariateData)]
            modelVar = modelVar[modelVar %in% colnames(covariateData)]
            modelString = paste0("coxph(patSurv ~ ", paste0(modelVar, collapse = " + "),  ", data=covariateData)")
            myModel = eval(parse(text=modelString))
            thisCoxDF = makeCoxDF(myModel, modelVar=modelVar, dfVar=myVar, returnGlobalStats = TRUE)

            # if ("years_to_birth" %in% colnames(covariateData)) {
            #     if ("gender" %in% colnames(covariateData)) {
            #         myModel = coxph(patSurv ~ years_to_birth + gender + meanMethyl + methylScore, data = covariateData)
            #         thisCoxDF = makeCoxDF(myModel, modelVar, dfVar, returnGlobalStats = TRUE)
            #     } else {
            #         myModel = coxph(patSurv ~ years_to_birth + meanMethyl + methylScore, data = covariateData)
            #     }
            # 
            # } else {
            #     if ("gender" %in% colnames(covariateData)) {
            #         myModel = coxph(patSurv ~ gender + meanMethyl + methylScore, data = covariateData)
            #     } else {
            #         myModel = coxph(patSurv ~ meanMethyl + methylScore, data = covariateData)
            #     }
            # 
            # }
            # add cancerID and regionSet name columns
            thisCoxDF = cbind(cancerID=rep(cancerID[i], nrow(thisCoxDF)), 
                              rsName=rep(abbrevName[j], nrow(thisCoxDF)), 
                              thisCoxDF)
            coxResults = rbind(coxResults, thisCoxDF)    
            

        sink(file = ffPlot(paste0(plotSubdir, "coxphModel", abbrevName[j], "_", cancerID[i], ".txt")))
        print(myModel)
        print(summary(myModel))
        sink()
        
        # test that model meets proportional hazards assumption
        # phTest = cox.zph(myModel)
        
        # cox model information
        # i is for cancer type, j is for region set, k is for Cox model variable 
        # myVar = c("methylScore", "age", "sex", "genomeMethyl", "globalStats")
        # modelVar = c("methylScore", "years_to_birth","gender", "meanMethyl", "GLOBAL")
        
       
        

        spearCorDF[2*(i-1)+j, "coxPVal"] = summary(myModel)$coefficients["methylScore", "Pr(>|z|)"]
        spearCorDF[2*(i-1)+j, "coxExp"] = summary(myModel)$coefficients["methylScore", "exp(coef)"]
        spearCorDF[2*(i-1)+j, "coxUpper"] = summary(myModel)$conf.int["methylScore", "upper .95"]
        spearCorDF[2*(i-1)+j, "coxLower"] = summary(myModel)$conf.int["methylScore", "lower .95"]
        # information about tests of cox model assumptions
        spearCorDF[2*(i-1)+j, "coxLower"] = summary(myModel)$conf.int["methylScore", "lower .95"]
        

            })
        # kaplan meier
        # create groups
        covariateData$methylGroup = ecdf(covariateData$methylScore)(covariateData$methylScore)
        covariateData$methylGroup[covariateData$methylGroup < 0.3333] = 0
        covariateData$methylGroup[covariateData$methylGroup > 0.6666] = 2
        covariateData$methylGroup[(covariateData$methylGroup >= 0.3333) & (covariateData$methylGroup <= 0.6666)] = NA

        kmFit = survfit(patSurv ~ methylGroup, data=covariateData)
        kmPlot = ggsurvplot(kmFit, risk.table = TRUE, conf.int = TRUE)

        # Error in data.frame(..., check.names = FALSE) :
        #     arguments imply differing number of rows: 165, 0, 330


        # cumcensor = TRUE, cumevents = TRUE,
        # a$plot = a$plot +
        # scale_color_discrete(name="Strata", labels=c("Low methylation score",
        #                                                    "High methylation score"), breaks=c("red", "blue"))
        #theme(legend.text = element_text(c("Low methylation score", "High methylation score")))
        pdf(file = ffPlot(paste0(plotSubdir, "kmPlot", abbrevName[j], "_", cancerID[i], ".pdf")))
        kmPlot
        dev.off()
        
        # also want to know how many covered CpGs are shared by EZH2 and SUZ12
        # only need to do once
        if (getOverlap) {

            pOL = percentCOverlap(mCoord = signalCoord, 
                                  GRList = myTopRS) 
            getOverlap = FALSE
        }
    }

}
coxResults = coxResults[-1, ] # blank first row

# spearCorDF = read.csv(ffSheets("EZH2_TCGA_cancer_stage.csv"))
spearCorDF$holmCoxPVal = p.adjust(p = spearCorDF$coxPVal, method = "holm")
spearCorDF$holmSpearmanPVal = p.adjust(p = spearCorDF$spearmanPVal, method = "holm")

write.csv(spearCorDF, file = ffSheets("polycomb_TCGA_cancer_stage.csv"), quote = FALSE, 
          row.names = FALSE, col.names = TRUE)
write.csv(coxResults, file = ffSheets("polycomb_TCGA_cancer_stage_coxPH.csv"), quote = FALSE, 
          row.names = FALSE, col.names = TRUE)

######################################################################################
# figures

# plotting confidence intervals for hazard ratios
# forest plot
# CompHome()
# coxResults=read.csv(file = "./processed/COCOA_paper/analysis/sheets/polycomb_TCGA_cancer_stage_coxPH.csv",
#          stringsAsFactors = FALSE)
coxResults$cancerID = factor(coxResults$cancerID, 
                             levels = rev(unique(coxResults$cancerID)))
fPlot <- ggplot(data=filter(coxResults, variableName=="methylScore"), 
                aes(x=cancerID, y=coxHRMean, ymin=coxHRLower, ymax=coxHRUpper)) +
    geom_pointrange() + 
    geom_hline(yintercept=1, lty=2) +
    coord_flip() +
    xlab("Cancer type") + ylab("Hazard ratio (95% CI)") +
    theme_classic() + scale_y_log10() 

fPlot
ggsave(filename = ffPlot(paste0(plotSubdir, "coxHazardRatios.svg")), 
       plot = fPlot, device = "svg")



# #################################################################################
# # are the same regions variable/associated with survival in each cancer?
# # are the same regions associated with cancer stage that are associated with survival?
# 
# 
# # get set of regions that is covered in each cancer
# 
# # get set of CpGs that are covered in all cancers
# 
# 
# # get association association with survival for each CpG, correcting for 
# # age, average genome-wide methylation, and sex
# cgNames = row.names(genomicSignal)
# covariateData = cbind(covariateData, t(genomicSignal))
# 
# try({
#     if ("years_to_birth" %in% colnames(covariateData)) {
#         if ("gender" %in% colnames(covariateData)) {
#             f1 <- as.formula(paste("patSurv ~ ",
#                                    paste(c("years_to_birth", "gender", "meanMethyl", cgNames), collapse= "+")))
#             
#             myModel = coxph(f1, covariateData)
#             myModel = coxph(patSurv ~ years_to_birth + gender + meanMethyl, 
#                             data = covariateData)
#         } else {
#             myModel = coxph(patSurv ~ years_to_birth + meanMethyl + methylScore, data = covariateData)
#         }
#         
#     } else {
#         if ("gender" %in% colnames(covariateData)) {
#             myModel = coxph(patSurv ~ gender + meanMethyl + methylScore, data = covariateData)
#         } else {
#             myModel = coxph(patSurv ~ meanMethyl + methylScore, data = covariateData)
#         }
#         
#     }
# })
# 
# coef(myModel)
# 
# 
# 
# makeQuantileGroups <- function(dataVec, nGroups=3) {
#     groups = rep(1, length(dataVec))
#     
#     for (i in 1:(nGroups-1)) {
#         groups[dataVec > quantile(x = dataVec, i/nGroups)] = i + 1
#     }
#     
#     # returns a vector of group membership
#     return(groups)
# }
# 
# 
# 




