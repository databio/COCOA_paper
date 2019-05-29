# does COCOA work with correlation values (instead of loadings)
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
# library(fastICA)

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients

# DNA methylation data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))

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

########################
# function to combine p values for features that overlap with a region set

combinePValRS <- function(dataDT, 
                          regionGR,
                          signalCols = colnames(dataDT)[!(colnames(dataDT) %in% c("chr", "start", "end"))]) { 
    total_region_number <- length(regionGR)
    mean_region_size <- round(mean(width(regionGR)), 1)
    dataGR <- COCOA:::BSdtToGRanges(list(dataDT))[[1]]
    OL <- findOverlaps(query = regionGR, subject = dataGR)
    if (length(OL) == 0) {
        return(NULL)
    }
    olCpG <- subjectHits(OL)
    # nonOLCpG <- (seq_len(nrow(dataDT)))[-olCpG]
    cytosine_coverage <- length(unique(olCpG))
    region_coverage <- length(unique(queryHits(OL)))
    combineP <- function(pval) {
        return(pchisq((sum(log(pval))*-2), df=length(pval)*2, lower.tail=F))
    }
    pVals <- vapply(X = signalCols, 
                    FUN = function(x) combineP(as.numeric(as.matrix(dataDT[olCpG, x, with = FALSE]))), FUN.VALUE = 1)
    wRes <- data.frame(t(pVals), cytosine_coverage, region_coverage, 
                       total_region_number, mean_region_size)
    return(wRes)
}


featurePCCorPVal = apply(X = pcScores[, PCsToAnnotate], MARGIN = 2, function(y) apply(X = centeredMData, 1, FUN = function(x) cor.test(x = x, y)$p.value))

hist(featurePC)
hist(-1* log(x = featurePCCorPVal[,"PC1"], base = 10))
hist(-1* log(x = featurePCCorPVal[,"PC2"], base = 10))
hist(-1* log(x = featurePCCorPVal[,"PC3"], base = 10))

##################################

featurePCCorPVal = featurePCCorPVal * nrow(featurePCCorPVal)
# maximum p val can be 1
featurePCCorPVal[featurePCCorPVal > 1] = 1

scoringMetric = "regionMean"

featurePCCorPVal = cbind(signalCoord, featurePCCorPVal)

simpleCache("combinePValScores", {
    rsScore = rbindlist(lapply(X = GRList[1:1000], FUN = function(x) combinePValRS(dataDT = featurePCCorPVal, regionGR = x)))
    rsScore$rsName = rsName
    rsScore$rsDescription= rsDescription
    rsScore
}, recreate=overwriteRSScoreResultsCaches)

View(rsScore[order(rsScore$PC1, decreasing = TRUE), ])
View(rsScore[order(rsScore$PC2, decreasing = TRUE), ])
View(rsScore[order(rsScore$PC4, decreasing = TRUE), ])
hist(rsScore$PC1)
hist(-1* log(x = rsScore$PC1, base = 10))



##############################################################################
# visualization pipeline






