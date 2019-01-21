# library(projectInit)

# project.init(codeRoot = paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/"), dataDir = paste0(Sys.getenv("PROCESSED"), "brca_PCA/"))
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
# library(fastICA)

# 
setwd(paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/plots/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)


# DNA methylation data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))

#############################################################################
# script specific IDs

plotSubdir = "02_5-brca_COCOA_cor/"
dataID = "657" # 657 patients with both ER and PGR info in metadata, 692 total
allMPCAString = "allMPCA_657" #  "allMPCA_657"
top10MPCAString = "top10MPCA_657"
rsScoreCacheName = paste0("rsScore_Cor_", dataID)
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


simpleCache(allMPCAString, assignToVariable = "allMPCA")
loadingMat = allMPCA$rotation
pcScores = allMPCA$x
signalCoord = brcaMList$coordinates


#############################################################################

cpgMeans = rowMeans(filteredMData)
# centering before calculating correlation
centeredMData = apply(X = filteredMData, MARGIN = 2, function(x) x - cpgMeans)
if (!all(colnames(filteredMData) == row.names(pcScores))) {
    stop("samples are not in the correct order in data structures.")
}

PCsToAnnotate = paste0("PC", 1:10)
# create feature correlation matrix with PCs (rows: features/CpGs, columns:PCs)
# how much do features correlate with each PC?
featurePCCor = apply(X = pcScores[, PCsToAnnotate], MARGIN = 2, function(y) apply(X = centeredMData, 1, FUN = function(x) cor(x = x, y)))

corLoadRatio = loadingMat[, PCsToAnnotate] / featurePCCor 
hist(corLoadRatio[, "PC10"])
# corLoadRatio = log(abs(loadingMat[, PCsToAnnotate]))/ featurePCCor 
# hist(corLoadRatio[(corLoadRatio[, "PC10"] <100) & (corLoadRatio[, "PC10"] > -100), "PC10"])

#############################################################################
# use featurePCCor instead of loadingMat in COCOA pipeline

scoringMetric = "regionMean"

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

View(rsScore[order(rsScore$PC1, decreasing = TRUE), ])
View(rsScore[order(rsScore$PC2, decreasing = TRUE), ])
View(rsScore[order(rsScore$PC4, decreasing = TRUE), ])
hist(rsScore$PC1)

##############################################################################
# visualization pipeline






