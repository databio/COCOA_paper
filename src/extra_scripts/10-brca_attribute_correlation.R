# correlate a sample attribute (like survival) with DNA methylation the run COCOA
# also try principal component regression

library(pls)
library(MultiAssayExperiment)

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))


#############################################################################
# script specific IDs

plotSubdir = "10-brca_attr_cor/"
dataID = "215" # 657 patients with both ER and PGR info in metadata, 692 total
rsScoreCacheName = paste0("rsScore_Surv_Cor_", dataID)
rsScoreCacheName2 = paste0("rsScore_Surv_PCR_", dataID)

overwriteRSScoreResultsCaches = TRUE

###############################################################################

simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList")

#############
#restrict patients included in this analysis
patientMetadata = patientMetadata[patientMetadata$subject_ID %in% 
                                      colnames(brcaMList[["methylProp"]]), ]
# getting rid of columns with no information on days of survival
patientMetadata = patientMetadata[(patientMetadata$vital_status == "alive" & !is.na(patientMetadata$days_to_last_follow_up)) | 
                                      (patientMetadata$vital_status == "dead" & !is.na(patientMetadata$days_to_death)), ]

# patientMetadata should have already screened out patients without ER/PGR status
# resulting in 657 patients
hasER_PGR_IDs = patientMetadata[, subject_ID]
filteredMData = brcaMList[["methylProp"]][, 
                                          colnames(brcaMList[["methylProp"]]) %in% hasER_PGR_IDs] 


###########################################################
# reading in the region sets
# load LOLA database
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))

#################################################################
# correlate DNA methylation with survival

days = seq(365, 365* 10, by = 0.5)
median(patientMetadata[patientMetadata$vital_status == "alive",]$days_to_last_follow_up, na.rm = TRUE)
hist(patientMetadata[patientMetadata$vital_status == "dead",]$days_to_death)
median(patientMetadata[patientMetadata$vital_status == "dead",]$days_to_death, na.rm = TRUE)
brcaSurv = patientMetadata$vital_status
sum(patientMetadata[patientMetadata$vital_status == "dead",]$days_to_death < 730)


nAlive = nrow(patientMetadata[patientMetadata$vital_status == "alive",])
totalDead = nrow(patientMetadata[patientMetadata$vital_status == "dead",])
numberDead = (ecdf(x = patientMetadata[patientMetadata$vital_status == "dead",]$days_to_death)(days) * totalDead)
numberKnownAlive = (1- ecdf(x = patientMetadata[patientMetadata$vital_status == "alive",]$days_to_last_follow_up)(days-1)) * nAlive
plot(days, 
     ecdf(x = patientMetadata[patientMetadata$vital_status == "dead",]$days_to_death )(days) * totalDead,  col = "red")
lines(days, 
     (1- ecdf(x = patientMetadata[patientMetadata$vital_status == "alive",]$days_to_last_follow_up)(days)) * nAlive, col="green")
# ratio of those known to be alive to those known to be dead
ratio = (numberKnownAlive + (totalDead - numberDead)) /
    numberDead

plot(days, ratio)

daysCutoff = days[ratio <= 4][1]
ratio[ratio <= 4][1]

# patients who died before daysCutoff are counted as dead, patients who are known
# to have lived past then are counted as alive
patientMetadata$keepCol = rep(0, nrow(patientMetadata))
# keep rows labeled with 1
patientMetadata$keepCol[patientMetadata$vital_status == "alive" & (patientMetadata$days_to_last_follow_up >= daysCutoff)] = 1
# some of these are known to be alive at daysCutoff and some are known to be dead, divide later
patientMetadata$keepCol[patientMetadata$vital_status == "dead"] = 1
patientMetadata = patientMetadata[patientMetadata$keepCol == 1, ]
patientMetadata$surv_at_daysCutoff = rep(0, nrow(patientMetadata))
patientMetadata$surv_at_daysCutoff[patientMetadata$vital_status == "alive"] = 1
patientMetadata$surv_at_daysCutoff[(patientMetadata$vital_status == "dead") & (patientMetadata$days_to_death > daysCutoff)] = 1
table(patientMetadata$surv_at_daysCutoff)
filteredMData = filteredMData[, as.character(patientMetadata$subject_ID)]


#### convert DNA methylation matrix to correlation matrix
# calculate correlation
featurePCCor = createCorFeatureMat(dataMat = filteredMData, 
                                   featureMat = as.matrix(patientMetadata$surv_at_daysCutoff), 
                                   centerDataMat=TRUE, centerFeatureMat=TRUE)
colnames(featurePCCor) <- "surv_at_daysCutoff"

#############################################################################
# run COCOA analysis

scoringMetric = "regionMean"
signalCoord = brcaMList$coordinates
PCsToAnnotate = "surv_at_daysCutoff"

simpleCache(rsScoreCacheName, {
    rsScore = runCOCOA(loadingMat=abs(featurePCCor), 
                       signalCoord = signalCoord, 
                       GRList, 
                       PCsToAnnotate = PCsToAnnotate, 
                       scoringMetric=scoringMetric)
    rsScore$rsName = rsName
    rsScore$rsDescription= rsDescription
    rsScore
}, recreate=overwriteRSScoreResultsCaches)

View(rsScore[order(rsScore$surv_at_daysCutoff, decreasing = TRUE), ])
hist(rsScore$surv_at_daysCutoff)

write.csv(x = rsScore, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/rsScore_methylCor_", dataID, ".csv"),
          quote = FALSE, row.names = FALSE)

##################################################################################
# now doing principal component regression and comparing with simple correlation

nPCs = 150
survSubsetPCA = prcomp(x = t(filteredMData), center = TRUE, scale. = FALSE)
varExpl = (survSubsetPCA$sdev^2 / sum(survSubsetPCA$sdev^2))
plot(varExpl[1:10])
sum(varExpl > 0.01)
survSubX = survSubsetPCA$x
survSubX = as.data.frame(cbind(survSubX[, 1:nPCs], surv_at_daysCutoff = patientMetadata$surv_at_daysCutoff))
a= lm(formula = surv_at_daysCutoff ~ ., data = survSubX)
pVals = summary(a)$pvals

# don't include intercept
weightedCor = abs(survSubsetPCA$rotation[, 1:nPCs]) %*% abs(coefficients(a)[-1]) 
weightedCor2 = (survSubsetPCA$rotation[, 1:nPCs]) %*% (coefficients(a)[-1]) 
hist(weightedCor2)
colnames(weightedCor) = "surv_at_daysCutoff"

#############################################################################
# run COCOA analysis
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))

scoringMetric = "regionMean"
signalCoord = brcaMList$coordinates
PCsToAnnotate = "surv_at_daysCutoff"

simpleCache(rsScoreCacheName2, {
    rsScore = runCOCOA(loadingMat=weightedCor, 
                       signalCoord = signalCoord, 
                       GRList, 
                       PCsToAnnotate = PCsToAnnotate, 
                       scoringMetric=scoringMetric)
    rsScore$rsName = rsName
    rsScore$rsDescription= rsDescription
    rsScore
}, recreate=overwriteRSScoreResultsCaches)

View(rsScore[order(rsScore$surv_at_daysCutoff, decreasing = TRUE), ])
hist(rsScore$surv_at_daysCutoff)

write.csv(x = rsScore, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/rsScore_brcaMethylPCR_", dataID, ".csv"),
          quote = FALSE, row.names = FALSE)

plotRSConcentration(rsScores = rsScore, scoreColName = "surv_at_daysCutoff", pattern = "esr|eralpha|eraalpha")
plotRSConcentration(rsScores = rsScore, scoreColName = "surv_at_daysCutoff", pattern = "suz12|ezh2|h3k27me|h3k9me")
plotRSConcentration(rsScores = rsScore, scoreColName = "surv_at_daysCutoff", pattern = "h3k27ac")

cor.test(as.numeric(as.factor(patientMetadata$ER_status)), patientMetadata$surv_at_daysCutoff)


##################################################################################
# Principal component regression

########## process patient metadata 
# filter out patients that don't have HER2 info
row.names(patientMetadata) <- patientMetadata$subject_ID
patientMetadata = patientMetadata[colnames(brcaMList[[2]]), ]
brcaMSE = SummarizedExperiment(brcaMList[[2]][, as.character(patientMetadata$subject_ID)], colData=patientMetadata)
brcaMSE = brcaMSE[, brcaMSE$her2_status %in% c("Positive", "Negative", "Indeterminate", "Equivocal")]

# assign numerical values based on HER2 status, positive=1, negative=-1, indeterminate/equivocal = 0
her2_numeric = rep(-99, nrow(colData(brcaMSE)))
her2_numeric[brcaMSE$her2_status == "Negative"] <- -1
her2_numeric[(brcaMSE$her2_status == "Indeterminate") | (brcaMSE$her2_status == "Equivocal")] <- 0
her2_numeric[(brcaMSE$her2_status == "Positive")] <- 1
brcaMSE$her2_numeric = her2_numeric
brcaMSE = brcaMSE[, brcaMSE$her2_numeric != -99]
# as.numeric(as.factor(brcaMSE$ER_status))

###########

test = as.data.frame(cbind(her2_numeric = brcaMSE$her2_numeric, t(assays(brcaMSE)[[1]][1:8000,])))
testPCR = pcr(her2_numeric ~ ., data = test, validation = "none")
summary(testPCR)
validationplot(testPCR)
dim(scores(testPCR))
class(scores(testPCR))
scores(testPCR)[1:5, 1:5]
plot(scores(testPCR)[, c(2,3)])


##########################################################################
# partial least squares regression

test = as.data.frame(cbind(her2_numeric = brcaMSE$her2_numeric, t(assays(brcaMSE)[[1]])[, 1:500]))
testPLSR = plsr(her2_numeric ~ ., data = test, validation = "CV")
summary(testPLSR)
validationplot(testPLSR)
dim(scores(testPLSR))
class(scores(testPLSR))
scores(testPLSR)[1:5, 1:5]
plot(scores(testPLSR)[, c(2,3)])

# PCA first to reduce dimensions then PLS
nPCs = ncol(her2PCA$x)
testPLSR = plsr(her2_numeric ~ ., 
                data = as.data.frame(cbind(her2PCA$x[, 1:nPCs], her2_numeric = brcaMSE$her2_numeric)), 
                validation = "CV")
summary(testPLSR)
validationplot(testPLSR)
scores(testPLSR)[1:5, 1:5]
RMSEP(testPLSR)
dim(coefficients(testPLSR))
plot(brcaMSE$her2_numeric, predict(object = testPLSR)[, , 1])

##########################################################################
# my own implementation of PCR
her2PCA = prcomp(x = t(assay(brcaMSE, 1)), center = TRUE, scale. = FALSE)
varExpl = (her2PCA$sdev^2 / sum(her2PCA$sdev^2))
plot(varExpl[1:10])
nPCs = sum(varExpl > 0.0025)
# nPCs = ncol(her2PCA$x) - 10
her2subX = as.data.frame(cbind(her2PCA$x[, 1:nPCs], her2_numeric = brcaMSE$her2_numeric))
a = lm(formula = her2_numeric ~ ., data = her2subX)
b= lm(formula = her2_numeric ~ ., data = her2subX[her2subX$her2_numeric != 0,])
summary(a)
summary(b)
predHer2 = predict(a, newdata = her2subX)
predHer2 = predict(b, newdata = her2subX[her2subX$her2_numeric != 0,])
plot(her2subX$her2_numeric[her2subX$her2_numeric != 0], predHer2)
ggplot(data = as.data.frame(cbind(her2subX, predHer2)), mapping = aes(x = her2_numeric, y = predHer2)) + geom_boxplot(aes(group=her2_numeric))
ggplot(data = as.data.frame(cbind(her2subX[her2subX$her2_numeric != 0, ], predHer2)), mapping = aes(x = her2_numeric, y = predHer2)) + geom_boxplot(aes(group=her2_numeric))
wilcox.test(x = predHer2[her2subX[her2subX$her2_numeric != 0, ]$her2_numeric == -1], 
            y= predHer2[her2subX[her2subX$her2_numeric != 0, ]$her2_numeric == 1])
# screen based on p values
pVals = summary(b)$coefficients[, "Pr(>|t|)"]
sum(pVals < (0.05/nPCs))
# get coef but don't include intercept
chosenPCcoef = coefficients(b)[-1][pVals[-1] < (0.05/nPCs)]
chosenPCs = names(chosenPCcoef)

#weightedLoad = abs(her2PCA$rotation[, chosenPCs]) %*% abs(chosenPCcoef) 
weightedLoad = abs(her2PCA$rotation[, chosenPCs] %*% (chosenPCcoef)) 
hist(weightedLoad)
colnames(weightedLoad) = "her2_weights"

#########
# visualize
ggplot(data = her2subX, mapping = aes(x = PC1, y = PC24)) + geom_point(aes(col = her2_numeric))
loadMed = apply(X = her2PCA$rotation, MARGIN = 2, median)
plot(loadMed)
loadRange = apply(X = her2PCA$rotation, MARGIN = 2, range)


######## 
# run COCOA on PCR-weighted loadings
# run COCOA analysis
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))

scoringMetric = "regionMean"
signalCoord = brcaMList$coordinates
PCsToAnnotate = "her2_weights"

simpleCache(paste0("her2PCRrsScores", "_", nPCs), {
    rsScore = runCOCOA(loadingMat=weightedLoad, 
                       signalCoord = signalCoord, 
                       GRList, 
                       PCsToAnnotate = PCsToAnnotate, 
                       scoringMetric=scoringMetric)
    rsScore$rsName = rsName
    rsScore$rsDescription= rsDescription
    rsScore
}, recreate=overwriteRSScoreResultsCaches)

View(rsScore[order(rsScore$her2_weights, decreasing = TRUE), ])
hist(rsScore$her2_weights)

write.csv(x = rsScore, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/rsScore_brcaMethylPCR_", dataID, ".csv"),
          quote = FALSE, row.names = FALSE)

plotRSConcentration(rsScores = rsScore, scoreColName = "her2_weights", pattern = "esr|eralpha|eraalpha")
plotRSConcentration(rsScores = rsScore, scoreColName = "her2_weights", pattern = "suz12|ezh2|h3k27me|h3k9me")
plotRSConcentration(rsScores = rsScore, scoreColName = "her2_weights", pattern = "h3k27ac")
plotRSConcentration(rsScores = rsScore, scoreColName = "her2_weights", pattern = "creb")

cor.test(as.numeric(as.factor(patientMetadata$ER_status)), patientMetadata$her2_weights)


#########
# example PCR code

data(iris)
library(pls)
a = pcr(Sepal.Length ~ ., data = iris, validation = "CV")
coef(a)
validationplot(a)
a$ncomp
dim(fitted(a))
plot(iris$Sepal.Length, fitted(a)[, 1, 1])
residuals(a)
dim(residuals(a))
plot(iris$Sepal.Length, residuals(a)[, 1, 1])
# https://www.rdocumentation.org/packages/pls/versions/2.7-1/topics/coef.mvr
