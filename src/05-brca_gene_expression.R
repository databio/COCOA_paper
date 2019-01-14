# showing that COCOA can integrate different data types
# first do PCA of gene expression, then correlate with DNA methylation
# do COCOA of DNA methylation


library("curatedTCGAData")

curatedTCGAData(diseaseCode = "BRCA", assays = "*", dry.run = TRUE)

rna = curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", 
                      dry.run = FALSE)
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


# correlate DNA methylation with gene expression