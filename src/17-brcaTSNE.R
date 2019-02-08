# does COCOA work with correlation and TSNE?

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(Rtsne)
# library(fastICA)

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)

# DNA methylation data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))

#############################################################################
# script specific IDs

plotSubdir = "17-brcaTSNE/"
dataID = "657" # 657 patients with both ER and PGR info in metadata, 692 total
allMPCAString = "allMPCA_657" #  "allMPCA_657"
top10MPCAString = "top10MPCA_657"
rsScoreCacheName = paste0("rsScore_TSNE_Cor_", dataID)
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
methVar = apply(filteredMData, 1, var)
ninetyQuant = quantile(methVar, probs = 0.9)
# top 10% variable CpGs
topFilteredMData = filteredMData[methVar >=ninetyQuant, ]
brcaTSNE = Rtsne(t(topFilteredMData))

#### convert DNA methylation matrix to correlation matrix
# calculate correlation
featurePCCor = createCorFeatureMat(dataMat = filteredMData, 
                                   featureMat = as.matrix(brcaTSNE$Y), 
                                   centerDataMat=TRUE, centerFeatureMat=TRUE)
colnames(featurePCCor) <- c("tSNE_Dim1", "tSNE_Dim2")


###########################################################
# reading in the region sets
# load LOLA database
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))

#############################################################################
# use featurePCCor instead of loadingMat in COCOA pipeline

scoringMetric = "regionMean"
PCsToAnnotate = c("tSNE_Dim1", "tSNE_Dim2")

simpleCache(rsScoreCacheName, {
    rsScore = runCOCOA(loadingMat=abs(featurePCCor), 
                       signalCoord = brcaMList$coordinates, 
                       GRList, 
                       PCsToAnnotate = PCsToAnnotate, 
                       scoringMetric=scoringMetric)
    rsScore$rsName = rsName
    rsScore$rsDescription= rsDescription
    rsScore
}, recreate=overwriteRSScoreResultsCaches)

View(rsScore[order(rsScore[, PCsToAnnotate[1]], decreasing = TRUE), ])
View(rsScore[order(rsScore[, PCsToAnnotate[2]], decreasing = TRUE), ])
View(rsScore[order(rsScore$PC4, decreasing = TRUE), ])
hist(rsScore$PC1)

##############################################################################
# visualization pipeline

plotRSConcentration(rsScores = rsScore, scoreColName = "tSNE_Dim2", pattern = "h3k4me1")

