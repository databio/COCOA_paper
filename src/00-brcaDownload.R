# download BRCA metadata

library("TCGAbiolinks")

library(data.table)
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
# http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html
# indexed clinical data is most up to date (data from GDCquery_clinic)
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

write.table(x = clinical, file = ffProjCode("metadata/brca_clinical_metadata.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 

# more in depth info is in the xml files
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical", file.type = "xml")

GDCdownload(query)
# available clinical info is given in package vignette
clinical1 <- GDCprepare_clinic(query, clinical.info = "patient")
clinical2 <- GDCprepare_clinic(query, clinical.info = "drug")
clinical3 <- GDCprepare_clinic(query, clinical.info = "follow_up")
clinical4 <- GDCprepare_clinic(query, clinical.info = "radiation")
clinical5 <- GDCprepare_clinic(query, clinical.info = "stage_event")
clinical6 <- GDCprepare_clinic(query, clinical.info = "new_tumor_event")
# there are also "admin" files
# write.csv(clinical1, file = "Other Clinical Data BRCA.csv")

# process BRCA metadata and connect with methylation data
bAnno = clinical1

# pulling out relevant info
bAnno$bcr_patient_barcode
bAnno$breast_carcinoma_progesterone_receptor_status
bAnno$breast_carcinoma_estrogen_receptor_status
bAnnoDT = data.table(subject_ID=bAnno$bcr_patient_barcode,
                     ER_status = bAnno$breast_carcinoma_estrogen_receptor_status,
                     ER_percent = bAnno$er_level_cell_percentage_category,
                     PGR_status = bAnno$breast_carcinoma_progesterone_receptor_status,
                     PGR_percent = bAnno$progesterone_receptor_level_cell_percent_category,
                     her2_status_IHC = bAnno$lab_proc_her2_neu_immunohistochemistry_receptor_status,
                     gender = bAnno$gender,
                     ethnicity = bAnno$ethnicity, 
                     race = bAnno$race_list,
                     menopause_status = bAnno$menopause_status,
                     age_diagnosis = bAnno$age_at_initial_pathologic_diagnosis,
                     her2_IHC_score = bAnno$her2_immunohistochemistry_level_result)
# previously her2_status_IHC listed as her2_status but changed for increased clarity

# get rid of any rows that are complete duplicates
bAnnoDT = unique(bAnnoDT)

write.csv(x = bAnnoDT, file = ffProjCode("metadata/brca_metadata.csv"), row.names = FALSE)

write.csv(x = clinical2, file = ffProjCode("metadata/brca_drug.csv"), row.names = FALSE)
