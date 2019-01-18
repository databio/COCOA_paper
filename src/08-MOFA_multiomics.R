# MOFA method
# public CLL data >>>>> hg19 <<<<<<
# this script relies on MOFA CLL vignette for the MOFA part of the analysis
library(ExperimentHub)
library("SummarizedExperiment")
library(MOFAtools)
library("FDb.InfiniumMethylation.hg19")


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
nrow(methCoordDT)
nrow(methData)

#ma

data("CLL_data")
MOFAobject <- createMOFAobject(CLL_data)

library(MultiAssayExperiment)


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

MOFAfactors <- getFactors(
    MOFAobject, 
    factors = c(1,2),
    as.data.frame = FALSE
)
head(MOFAfactors)
plot(MOFAfactors)

# get coordinates for CpG microarray probes (hg19)



