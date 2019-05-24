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

library(data.table)

# process BRCA metadata and connect with methylation data
brcaAnnoPath = "C:/Users/John L/Desktop/AML/brca/Other_Clinical_Data_BRCA.csv"
bAnno = fread(file = brcaAnnoPath)

# pulling out relevant info
bAnno$bcr_patient_barcode
bAnno$breast_carcinoma_progesterone_receptor_status
bAnno$breast_carcinoma_estrogen_receptor_status
bAnnoDT = data.table(subject_ID=bAnno$bcr_patient_barcode,
                     ER_status = bAnno$breast_carcinoma_estrogen_receptor_status,
                     PGR_status = bAnno$breast_carcinoma_progesterone_receptor_status,
                     her2_status = bAnno$lab_proc_her2_neu_immunohistochemistry_receptor_status,
                     gender = bAnno$gender,
                     ethnicity = bAnno$ethnicity, 
                     race = bAnno$race_list,
                     menopause_status = bAnno$menopause_status,
                     age_diagnosis = bAnno$age_at_initial_pathologic_diagnosis)
fileName = paste0(Sys.getenv("CODE"), 
                  "PCARegionAnalysis/metadata/brca_metadata.csv")
write.csv(x = bAnnoDT, file = fileName, row.names = FALSE)
