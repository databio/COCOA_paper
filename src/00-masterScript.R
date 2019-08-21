# This script runs the other R script files in order to produce the analysis

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

##############################################################################
# BRCA DNA methylation analysis 
source(ffProjCode("19-permTestBRCADNAmPCA.R"))
source(ffProjCode("03-brca_DNAm_COCOA_visualization")) # self-contained
source(ffProjCode("04-brca_DNAm_figures.R"))

##############################################################################
# BRCA ATAC analysis
source(ffProjCode("19-permTestBRCA_ATAC_PCA.R"))
source(ffProjCode("06_2-brca_ATAC_vis")) # self-contained (doesn't run main visualization pipeline)

##############################################################################
# MOFA analysis
source(ffProjCode("19-permTestMOFACLLDNAm.R"))
source(ffProjCode("09-MOFA_visualization.R")) # self-contained (loads data)
source(ffProjCode("09-MOFA_figures.R"))

##############################################################################
# KIRC analysis
source(ffProjCode("22-otherCancerAttrCorrelation.R")) # KIRC analysis
# KIRC, haven't run comprehensive visualization
source(ffProjCode("23-survivalAnalysis")) # KIRC, self-cont?
source(ffProjCode("22-KIRC_figures.R"))


###############################################################################
# creating additional figures


# source(ffProjCode(""))
