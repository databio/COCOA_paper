# correlate a sample attribute (like survival) with DNA methylation the run COCOA


source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)


#############################################################################


# try with her2 status
table(patientMetadata$her2_status)

#### convert DNA methylation matrix to correlation matrix




# calculate correlation with each latent factor from MOFA
featurePCCor = createCorFeatureMat(dataMat = methData, 
                                   featureMat = latentFactors, 
                                   center=TRUE)

##################################################################
# load hg38 region set database



##################################################################
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