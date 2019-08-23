# showing that COCOA can integrate different data types
# first do PCA of gene expression, then correlate with DNA methylation
# do COCOA of DNA methylation

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library("curatedTCGAData")

curatedTCGAData(diseaseCode = "BRCA", assays = "*", dry.run = TRUE)

brca = curatedTCGAData(diseaseCode = "BRCA", assays = c("Methylation*"), dry.run = FALSE)
testM = assays(brca)[[1]]
testM = as.matrix(testM)
testM[221:226, 934:937]
which(is.na(testM))


rna = curatedTCGAData(diseaseCode = "BRCA", assays = c("RNASeq2GeneNorm"), 
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

# pca of gene expression
rnaMat = assays(rna)[[1]]
rpca = prcomp(t(rnaMat))
plot(rpca$x[,"PC1"], rpca$x[, "PC4"])
wideFormat(DataFrame(rpca$x), 
           colDataCols = "patient.breast_carcinoma_estrogen_receptor_status")
erStatus = colData(rna)$patient.breast_carcinoma_estrogen_receptor_status
rpcaX = cbind(rpca$x, erStatus) 
colorClusterPlots()
nrow(rpca$x)
length(erStatus)
View(colData(rna))
length(unique(colnames(rnaMat)))
setdiff(colnames(rnaMat), colData(rna)$patientID)
# correlate DNA methylation with gene expression

