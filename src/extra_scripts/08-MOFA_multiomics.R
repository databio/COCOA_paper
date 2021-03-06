# MOFA method
# public CLL data >>>>> hg19 <<<<<<
# this script relies on MOFA CLL vignette for the MOFA part of the analysis
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(ExperimentHub)
library("SummarizedExperiment")
library("FDb.InfiniumMethylation.hg19")
library(MultiAssayExperiment)
# library(COCOA)
devtools::load_all(ffCode("COCOA/"))
library(MOFAdata)
library(MOFA)

#############################################################################
# script specific IDs

plotSubdir = "08-MOFA_multiomics/"
dataID = "CLL196" # 657 patients with both ER and PGR info in metadata, 692 total
variationType = "cov" # covariance instead of correlation
rsScoreCacheName = paste0("rsScore_", dataID, "MOFA", "_", variationType)
overwriteRSScoreResultsCaches = TRUE


if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

###############################################################################


cllMethyl = prepareCLLMethyl()
methData = cllMethyl$methylProp
methCoord = cllMethyl$methylCoord

# multiomics data for MOFA
data("CLL_data") 
# MOFAobject <- createMOFAobject(CLL_data)
# 
data("CLL_covariates")
head(CLL_covariates)

# Create MultiAssayExperiment object
mae_CLL <- MultiAssayExperiment(
    experiments = CLL_data,
    colData = CLL_covariates
)



# Build the MOFA object
MOFAobject <- createMOFAobject(mae_CLL)

# Loading an existing trained model
filepath <- system.file("extdata", "CLL_model.hdf5",
                        package = "MOFAdata")

MOFAobject <- loadModel(filepath, MOFAobject)

r2 <- calculateVarianceExplained(MOFAobject)
r2$R2Total


MOFAweights <- getWeights(
    MOFAobject, 
    views = "Methylation", 
    factors = "all", 
    as.data.frame = TRUE
)
head(MOFAweights)


# correlate DNA methylation and latent factor scores 

latentFactors <- getFactors(
    MOFAobject,
    as.data.frame = FALSE
)
head(latentFactors)
plot(latentFactors)

# make sure samples are in same order in DNA methylation data from MOFA and 
# CLL data package
# there are less samples in methylation data
latentFactors = latentFactors[colnames(methData), ]
if (!all(colnames(methData) == row.names(latentFactors))) {
    stop("Sample order is not the same in objects")
}

apply(X = latentFactors, MARGIN = 2, FUN = function(x) any(is.na(x)))
simpleCache("cllMOFAFactors", {cllMOFAFactors = latentFactors})

##################################################################
# load hg19 region set database
loadGRList(genomeV = "hg19")

#####################################################################

loadMOFAData(methylMat = FALSE, signalCoord=FALSE, latentFactorMat = FALSE,
             cllMultiOmics=TRUE)

mutDF = cllMultiOmics$Mutation
# same order
mutDF = mutDF[colnames(mutDF), ]
IGHV = mutDF["IGHV", ]
IGHV[IGHV == 1] = "mutated"
IGHV[IGHV == 0] = "unmutated"
IGHV = as.factor(IGHV)

######################################################################
# analysis of gender and latent factors
gata3RS =  GRList[["wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak"]]
# GRList[["wgEncodeAwgTfbsHaibT47dEraaV0416102Bpa1hUniPk.narrowPeak"]]#
# runCOCOA(signal = )
sharedSamples = row.names(CLL_covariates)[row.names(CLL_covariates) %in% colnames(methData)]
IGHV = mutDF[, sharedSamples]["IGHV", ]
unmutSamples = sharedSamples[IGHV==0 & !is.na(IGHV)]
g3PCA = dimRedOnRS(regionSet = gata3RS, methylData = methData[, unmutSamples], mCoord = methCoord, drMethod = "pca")
g3X = g3PCA$x
g3X = cbind(g3X[unmutSamples, ], CLL_covariates[unmutSamples, ], t(mutDF[,unmutSamples]))
# g3X = cbind(g3X[sharedSamples, ], CLL_covariates[sharedSamples, ], t(mutDF[,sharedSamples]))
colorClusterPlots(g3X, plotCols = paste0("PC", c(1,2)), colorByCols = c("Gender", "IGHV"))
maleID = row.names(CLL_covariates[CLL_covariates$Gender=="m",]) 
femaleID = row.names(CLL_covariates[CLL_covariates$Gender=="f",])
whichPC = "PC1"
whichPC = "PC2"
wilcox.test(g3X[maleID[maleID %in% row.names(g3X)], whichPC], g3X[femaleID[femaleID %in% row.names(g3X)], whichPC])

whichLF = "LF7"
whichLF = "LF1"
whichLF = "LF9"
wilcox.test(latentFactors[maleID[maleID %in% row.names(latentFactors)], whichLF], latentFactors[femaleID[femaleID %in% row.names(latentFactors)], whichLF])

plotFactorScatters(MOFAobject, factors = c(1,7, 9), color_by = "Gender")

######################################################################
#### convert DNA methylation matrix to correlation matrix (or covariance matrix)


# calculate correlation with each latent factor from MOFA
simpleCache(paste0("inferredMethylWeightsMOFA", "_", variationType), {
    featurePCCor = createCorFeatureMat(dataMat = methData, 
                                       featureMat = latentFactors, 
                                       centerDataMat = TRUE, 
                                       centerFeatureMat = TRUE, testType = variationType)
    inferredMethylWeightsMOFA = featurePCCor
})



##################################################################
# run COCOA analysis

scoringMetric = "regionMean"
signalCoord = methCoord
# # only latent factors with no NA's
# PCsToAnnotate = paste0("LF", c(1:3, 5:7, 9))
# include latent factors with NA's
PCsToAnnotate = paste0("LF", 1:10)

simpleCache(rsScoreCacheName, {
    rsScore = runCOCOA(signal=abs(featurePCCor), 
                       signalCoord = signalCoord, 
                       GRList, 
                       signalCol = PCsToAnnotate, 
                       scoringMetric=scoringMetric,
                       signalCoordType = "singleBase")
    rsScore$rsName = rsName
    rsScore$rsDescription= rsDescription
    rsScore
}, recreate=overwriteRSScoreResultsCaches, assignToVariable = "rsScore")
rsScores = as.data.table(rsScore)

View(rsScore[order(rsScore$LF1, decreasing = TRUE), ])
View(rsScore[order(rsScore$LF2, decreasing = TRUE), ])
View(rsScore[order(rsScore$LF4, decreasing = TRUE), ])
View(rsScore[order(rsScore$LF7, decreasing = TRUE), ])
hist(rsScore$LF7)

write.csv(x = rsScore, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/rsScore_", "methyl_"variationType "_", dataID, ".csv"),
          quote = FALSE, row.names = FALSE)


####################################################################
# permutation based analysis
featurePCCor = inferredMethylWeightsMOFA
signalCoord = methCoord
# only latent factors with no NA's
PCsToAnnotate = paste0("LF", c(1:3, 5:7, 9))


nPerm = 10000
sampleSize = c(10, 100, 1000, 10000, 100000)
simpleCache(paste0("nullDistList", dataID), {
    nullDistList = multiNullDist(loadingMat = featurePCCor, PCsToAnnotate=PCsToAnnotate, nPerm = nPerm, sampleSize = sampleSize)
    nullDistList
    
}, assignToVariable = "nullDistList")

for (j in seq_along(PCsToAnnotate)) {
    pcOfInt = PCsToAnnotate[j]
    pdf(ffPlot(paste0(plotSubdir, "nullDist", pcOfInt, "_", nPerm, "Perm.pdf"))) 
    for (i in seq_along(nullDistList)) {
        hist(nullDistList[[i]][, pcOfInt], main = paste0("sampleSize=", sampleSize[i]))    
    }
    dev.off()
}


getPValRangeVec <- function(rsScore, nullDistList, pc, regionCoverage) {
    getPValRange(rsScore, nullDistList, pc, regionCoverage)
}

# convert to scores (regionMean) to p values
upperBoundPVal = list()
for (i in seq_along(PCsToAnnotate)) {
    upperBoundPVal[[i]] = getUpperBound(rsScore = as.numeric(as.matrix(rsScores[, PCsToAnnotate[i], with=FALSE])), 
                                        nullDistList = nullDistList, 
                                        sampleSize=sampleSize, 
                                        pc=PCsToAnnotate[i], 
                                        regionCoverage = rsScores$region_coverage)
}

upperBoundPVal = do.call(cbind, upperBoundPVal)
colnames(upperBoundPVal) <- PCsToAnnotate
upperBoundPVal = as.data.table(upperBoundPVal)

upperBoundPVal = cbind(upperBoundPVal, rsScores[, c("cytosine_coverage", "region_coverage", 
                                                    "total_region_number", "mean_region_size",
                                                    "rsName", "rsDescription")])
colnames(upperBoundPVal) <- paste0(colnames(upperBoundPVal), "_pval")
rsScores = cbind(rsScores, upperBoundPVal)
i=1
setorderv(rsScores, cols = c(paste0(PCsToAnnotate[i], "_pval"), PCsToAnnotate[i]), order = c(1L, -1L))
plot(as.numeric(as.matrix(rsScores[, PCsToAnnotate[i], with=FALSE])))
plot(-log(as.numeric(as.matrix(rsScores[, paste0(PCsToAnnotate[i], "_pval"), with=FALSE]))))
View(rsScores)

# pdf(ffPlot(paste0(plotSubdir, "rsScoresPValHist", nPerm, "PermBRCA657.pdf")))
# for (i in seq_along(PCsToAnnotate)) {
#     hist(-log(as.numeric(as.matrix(upperBoundPVal[, PCsToAnnotate[i], with=FALSE])) * 2246, base = 10), main = paste0("-log10 p value for ", PCsToAnnotate[i]), xlab = "-log10 P Val")
# }
# 
# dev.off()
