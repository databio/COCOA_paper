# RPPA protein data


source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))


#############################################################################
# reading in the region sets
# load LOLA database
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))

#################################################################

simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList")



# proteins of interest (we have a region set related to them):
# NFKB, GATA3, GATA6, ERALPHA_pS118, AR, CJUN_pS73, CMYC, , 
# BETACATENIN, PR, SMAD1, SMAD3, STAT3_pY705, STAT5ALPHA, ETS1, FOXM1, IRF1, 
# interesting but not necessarily clear what the "right" answer is: 
# CD31, BRAF, CYCLIND1, RICTOR, CD49B, CD20
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

sum(colnames(brcaMList[[2]]) %in% prot$subject_ID)
sharedSamples = colnames(brcaMList[[2]])[colnames(brcaMList[[2]]) %in% prot$subject_ID]

filteredMData = brcaMList[[2]][, sharedSamples]
filtProt = prot[sharedSamples, ]

hist(filtProt$NFKBP65_pS536)
hist(filtProt$GATA3)

#### convert DNA methylation matrix to correlation matrix
# calculate correlation
featurePCCor = createCorFeatureMat(dataMat = filteredMData, 
                                   featureMat = as.matrix(filtProt[, proteinOfInterest]), 
                                   centerDataMat=TRUE, centerFeatureMat=TRUE)
colnames(featurePCCor) <- proteinOfInterest
hist(featurePCCor[, "GATA3"])

#############################################################################
# run COCOA analysis

scoringMetric = "regionMean"
signalCoord = brcaMList$coordinates
PCsToAnnotate = proteinOfInterest

simpleCache("rppa_cor_rsScores", {
    rsScore = runCOCOA(loadingMat=abs(featurePCCor), 
                       signalCoord = signalCoord, 
                       GRList, 
                       PCsToAnnotate = PCsToAnnotate, 
                       scoringMetric=scoringMetric)
    rsScore$rsName = rsName
    rsScore$rsDescription= rsDescription
    rsScore
}, recreate=TRUE)

View(rsScore[order(rsScore$GATA3, decreasing = TRUE), ])
hist(rsScore$GATA3)
dataID = ""
write.csv(x = rsScore, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/rsScore_methylCor_RPPA", dataID, ".csv"),
          quote = FALSE, row.names = FALSE)

###########################
# visualization
plotRSConcentration(rsScores = rsScore, scoreColName = "NFKBP65_pS536", 
                    colsToSearch = c("rsName", "rsDescription"), pattern = "nfkb|nfkappa")
plotRSConcentration(rsScores = rsScore, scoreColName = "CMYC", 
                    colsToSearch = c("rsName", "rsDescription"), pattern = "myc")
plotRSConcentration(rsScores = rsScore, scoreColName = "STAT3_pY705", 
                    colsToSearch = c("rsName", "rsDescription"), pattern = "eralpha|eraa|esr1")
plotRSConcentration(rsScores = rsScore, scoreColName = "STAT3_pY705", 
                    colsToSearch = c("rsName", "rsDescription"), pattern = "stat3")
plotRSConcentration(rsScores = rsScore, scoreColName = "GATA3", 
                    colsToSearch = c("rsName", "rsDescription"), pattern = "gata3")
plotRSConcentration(rsScores = rsScore, scoreColName = "ERALPHA_pS118", 
                    colsToSearch = c("rsName", "rsDescription"), pattern = "eralpha|eraa|esr1")
plotRSConcentration(rsScores = rsScore, scoreColName = "SMAD3", 
                    colsToSearch = c("rsName", "rsDescription"), pattern = "eralpha|eraa|esr1")
plotRSConcentration(rsScores = rsScore, scoreColName = "SMAD3", 
                    colsToSearch = c("rsName", "rsDescription"), pattern = "smad3")
absMedCor = apply(abs(featurePCCor), MARGIN = 2, FUN = median)
absMaxCor = apply(abs(featurePCCor), MARGIN = 2, FUN = max)
hist(absMedCor)

##################################################################################
# try principal component regression instead
# my own implementation of PCR
filtProt
allMPCAString = "allMPCA_657"
simpleCache(allMPCAString, assignToVariable = "mPCA")
varExpl = (mPCA$sdev^2 / sum(mPCA$sdev^2))
plot(varExpl[1:10])
nPCs = sum(varExpl > 0.0025)
# nPCs = ncol(mPCA$x) - 10
sharedPCASamples = rownames(filtProt) %in% rownames(mPCA$x)
mPCASubX = as.data.frame(cbind(mPCA$x[rownames(filtProt)[sharedPCASamples], 1:nPCs], proteinLevel = filtProt$STAT5ALPHA[sharedPCASamples]))
a = lm(formula = proteinLevel ~ ., data = mPCASubX)
summary(a)
predHer2 = predict(a, newdata = mPCASubX)
plot(predHer2, mPCASubX$proteinLevel)
# screen based on p values
pVals = summary(a)$coefficients[, "Pr(>|t|)"]
sum(pVals < 0.05)
sum(pVals < (0.05/nPCs))
# get coef but don't include intercept
chosenPCcoef = coefficients(b)[-1][pVals[-1] < (0.05/nPCs)]
chosenPCs = names(chosenPCcoef)

#weightedLoad = abs(her2PCA$rotation[, chosenPCs]) %*% abs(chosenPCcoef) 
weightedLoad = abs(her2PCA$rotation[, chosenPCs] %*% (chosenPCcoef)) 
hist(weightedLoad)
colnames(weightedLoad) = "proteinLevel"