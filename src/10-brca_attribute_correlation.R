# correlate a sample attribute (like survival) with DNA methylation the run COCOA


source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))


#############################################################################
# script specific IDs

plotSubdir = "10-brca_attr_cor/"
dataID = "657" # 657 patients with both ER and PGR info in metadata, 692 total
allMPCAString = "allMPCA_657" #  "allMPCA_657"
top10MPCAString = "top10MPCA_657"
rsScoreCacheName = paste0("rsScore_Surv_Cor_", dataID)
overwriteRSScoreResultsCaches = TRUE

###############################################################################

simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList")
#restrict patients included in this analysis
patientMetadata = patientMetadata[patientMetadata$subject_ID %in% 
                                      colnames(brcaMList[["methylProp"]]), ]
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
numberKnownAlive = (1- ecdf(x = patientMetadata[patientMetadata$vital_status == "alive",]$days_to_last_follow_up)(days)) * nAlive
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
patientMetadata$keepCol[patientMetadata$vital_status == "dead" & (patientMetadata$days_to_death <= daysCutoff)] = 1
patientMetadata = patientMetadata[patientMetadata$keepCol == 1, ]
table(patientMetadata$vital_status)

#### convert DNA methylation matrix to correlation matrix
# calculate correlation
featurePCCor = createCorFeatureMat(dataMat = methData, 
                                   featureMat = as.matrix(brcaSurv), 
                                   centerDataMat=TRUE, centerFeatureMat=TRUE)

#############################################################################
# run COCOA analysis

scoringMetric = "regionMean"
signalCoord = methCoord
# only latent factors with no NA's
PCsToAnnotate = paste0("LF", c(1:3, 5:7, 9))

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

View(rsScore[order(rsScore$LF1, decreasing = TRUE), ])
View(rsScore[order(rsScore$LF2, decreasing = TRUE), ])
View(rsScore[order(rsScore$LF4, decreasing = TRUE), ])
View(rsScore[order(rsScore$LF7, decreasing = TRUE), ])
hist(rsScore$LF7)

write.csv(x = rsScore, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/rsScore_methylCor_", dataID, ".csv"),
          quote = FALSE, row.names = FALSE)


#