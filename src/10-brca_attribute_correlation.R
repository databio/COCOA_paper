# correlate a sample attribute (like survival) with DNA methylation the run COCOA
library(pls)

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