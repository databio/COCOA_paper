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

setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))

#############################################################################
# script specific IDs

plotSubdir = "08-MOFA_multiomics/"
dataID = "CLL196" # 657 patients with both ER and PGR info in metadata, 692 total
rsScoreCacheName = paste0("rsScore_Cor_", dataID)
overwriteRSScoreResultsCaches = TRUE

###############################################################################


# match probe names with genomic coordinates
ls('package:FDb.InfiniumMethylation.hg19')
m450kAnno = get450k()
length(m450kAnno)

eh = ExperimentHub()
query(eh, "CLLmethylation")

meth = eh[["EH1071"]]
# rows are cpgs, columns are samples
methData = assay(meth)
dataProbeNames = row.names(methData)
# methData[1:5, 1:5]
# get coordinates in same order as CLL data
methCoord = m450kAnno[dataProbeNames]
all(names(methCoord) == dataProbeNames)
methCoordDT = COCOA:::grToDt(methCoord)
# keep start coordinate as CpG site
methCoordDT = methCoordDT[, .(chr, start)]
methCoord = COCOA:::dtToGr(methCoordDT)
nrow(methCoordDT)
nrow(methData)

#ma

data("CLL_data")
MOFAobject <- createMOFAobject(CLL_data)




data("CLL_covariates")
head(CLL_covariates)

# Create MultiAssayExperiment object 
mae_CLL <- MultiAssayExperiment(
    experiments = CLL_data, 
    colData = CLL_covariates
)

# Build the MOFA object
MOFAobject <- createMOFAobject(mae_CLL)

plotTilesData(MOFAobject)

DataOptions <- getDefaultDataOptions()
DataOptions 

ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 25
ModelOptions


TrainOptions <- getDefaultTrainOptions()

# Automatically drop factors that explain less than 2% of variance in all omics
TrainOptions$DropFactorThreshold <- 0.02

TrainOptions$seed <- 2017

TrainOptions


MOFAobject <- prepareMOFA(
    MOFAobject, 
    DataOptions = DataOptions,
    ModelOptions = ModelOptions,
    TrainOptions = TrainOptions
)


MOFAobject <- runMOFA(MOFAobject)

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

#### convert DNA methylation matrix to correlation matrix

# @param dataMat columns of dataMat should be samples/patients, rows should be genomic signal
# (each row corresponds to one genomic coordinate/range)
# @param featureMat Rows should be samples, columns should be "features" 
# (whatever you want to get correlation with: eg PC scores),
# all columns in featureMat will be used (subset when passing to function
# in order to not use all columns)
# @param center logical object. Should rows in dataMat be centered based on
# their means? (subtracting row means from each row)
#
# returns a matrix where rows are the genomic signal (eg a CpG or region) and
# columns are the columns of featureMat
createCorFeatureMat = function(dataMat, featureMat, center=TRUE) {
    if (center) {
        cpgMeans = rowMeans(dataMat)
        # centering before calculating correlation
        dataMat = apply(X = dataMat, MARGIN = 2, function(x) x - cpgMeans)
        
    }

    
    # create feature correlation matrix with PCs (rows: features/CpGs, columns:PCs)
    # how much do features correlate with each PC?
    featurePCCor = apply(X = featureMat, MARGIN = 2, function(y) apply(X = dataMat, 1, FUN = function(x) cor(x = x, y)))
    return(featurePCCor)
    # corLoadRatio = loadingMat[, PCsToAnnotate] / featurePCCor 
    # hist(corLoadRatio[, "PC10"])
}


# calculate correlation with each latent factor from MOFA
featurePCCor = createCorFeatureMat(dataMat = methData, 
                                   featureMat = latentFactors, 
                                   center=TRUE)

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


