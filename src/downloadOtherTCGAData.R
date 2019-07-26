
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(curatedTCGAData)
library(TCGAutils)

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)

#############################################################################

cancerID = "KIRC"
#########

curatedTCGAData(diseaseCode = cancerID, assays = "*", dry.run = TRUE)

mData = curatedTCGAData(diseaseCode = cancerID, assays = c("*methyl450*"), dry.run = FALSE)
testM = assays(mData)[[1]]
testM = as.matrix(testM)
testM[221:226, 304:307]
hasNA = apply(testM, 2, function(x) sum(is.na(x)))
hist(hasNA)
cpgHasNA = apply(testM, 1, function(x) sum(is.na(x)))
table(cpgHasNA)


##### get CpG coordinates
# match probe names to coordinates
# screen out CpGs that have any NAs and the XY chromosomes
methylList = filtMethylMat(signalCoord=NULL, methylMat = testM) 



############
# get patient metadata, stored in colData(curatedTCGAData())
metaDataCols = getClinicalNames(cancerID)
allMeta = colData(mData) 
tcgaMetadata = allMeta[, metaDataCols]



############

rna = curatedTCGAData(diseaseCode = cancerID, assays = c("RNASeq2GeneNorm"), 
                      dry.run = FALSE)
dlbc <- curatedTCGAData("DLBC", assays = c("RNASeq2GeneNorm", "Methylation"), FALSE)
class(rna)
assays(rna)[[1]][1:5, 1:5]

# subset all data down to samples that were previously used in the DNA
# methylation analysis
test = rna[, 1:20, ]
head(assays(test)[[1]])
test$patientID
test@colData@listData
sampleMap(test)
testMD = colData(test)
testMD2 = as.data.frame(testMD)
head(testMD2)
grep(pattern = "estrogen_receptor_status", colnames(testMD), value=TRUE)
wideFormat(test, 
           colDataCols = "patient.breast_carcinoma_estrogen_receptor_status")