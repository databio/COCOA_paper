# This script runs the other R script files in order to produce the analysis

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

###############################################################################
# data preparation


###############################################################################
# COCOA analysis

source(ffProjCode("19-permTestBRCADNAmPCA.R"))
source(ffProjCode("19-permTestBRCA_ATAC_PCA.R"))
source(ffProjCode("19-permTestMOFACLLDNAm.R"))
source(ffProjCode("22-otherCancerAttrCorrelation.R")) # KIRC analysis



###############################################################################
# creating figures