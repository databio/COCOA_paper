# download BRCA metadata

library("TCGAbiolinks")

library(dplyr)
library(readr)
library(SummarizedExperiment)
tcgaDataPath = paste0(getOption("PROCESSED"), "/", "tcga")
if (!dir.exists(tcgaDataPath)) {
    dir.create(tcgaDataPath)
}
setwd(tcgaDataPath)

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))


###################################################

# allBiospecimen = GDCquery("TCGA-BRCA", data.category = "Biospecimen")
# allBiospecimenDF = getResults(allBiospecimen)
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

write.table(x = clinical, file = dirCode("metadata/brca_clinical_metadata.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 


##############################################################

