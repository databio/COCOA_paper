# MOFA method
# public CLL data >>>>> hg19 <<<<<<
# this script relies on MOFA CLL vignette for the MOFA part of the analysis
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(ExperimentHub)
library("SummarizedExperiment")
library(MOFAtools)
library("FDb.InfiniumMethylation.hg19")
library(MultiAssayExperiment)
library(COCOA)

setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))

#############################################################################
# script specific IDs

plotSubdir = "08-MOFA_multiomics/"
dataID = "CLL196" # 657 patients with both ER and PGR info in metadata, 692 total
rsScoreCacheName = paste0("rsScore_Cor_", dataID)
overwriteRSScoreResultsCaches = TRUE

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
# 
# plotTilesData(MOFAobject)
# 
# DataOptions <- getDefaultDataOptions()
# DataOptions 
# 
# ModelOptions <- getDefaultModelOptions(MOFAobject)
# ModelOptions$numFactors <- 25
# ModelOptions
# 
# 
# TrainOptions <- getDefaultTrainOptions()
# 
# # Automatically drop factors that explain less than 2% of variance in all omics
# TrainOptions$DropFactorThreshold <- 0.02
# 
# TrainOptions$seed <- 2017
# 
# TrainOptions
# 
# 
# MOFAobject <- prepareMOFA(
#     MOFAobject, 
#     DataOptions = DataOptions,
#     ModelOptions = ModelOptions,
#     TrainOptions = TrainOptions
# )
# 
# 
# MOFAobject <- runMOFA(MOFAobject)

# Loading an existing trained model
filepath <- system.file("extdata", "CLL_model.hdf5",
                        package = "MOFAtools")

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

#####################################################################\
# analysis of gender and latent factors
maleID = row.names(CLL_covariates[CLL_covariates$Gender=="m",]) 
femaleID = row.names(CLL_covariates[CLL_covariates$Gender=="f",])
whichLF = "LF7"
wilcox.test(latentFactors[maleID, whichLF], latentFactors[femaleID, whichLF])

plotFactorScatters(MOFAobject, factors = c(7, 9), color_by = "Gender")

######################################################################
#### convert DNA methylation matrix to correlation matrix

# calculate correlation with each latent factor from MOFA
simpleCache("inferredMethylWeightsMOFA", {
    featurePCCor = createCorFeatureMat(dataMat = methData, 
                                       featureMat = latentFactors, 
                                       centerDataMat = TRUE, 
                                       centerFeatureMat = TRUE)
    inferredMethylWeightsMOFA = featurePCCor
})


##################################################################
# load hg19 region set database

source(paste0(Sys.getenv("CODE"), "aml_e3999/src/load_process_regionDB.R"))

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


