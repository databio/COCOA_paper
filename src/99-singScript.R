# script to run on singularity on cluster

# library for singularity
.libPaths(c("/scratch/jtl2hk/3.6", .libPaths()))

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

#############################################################################
# source(ffProjCode("19-permTestBRCADNAmPCA.R"))
# source(ffProjCode("03-brca_DNAm_COCOA_visualization")) # self-contained
# source(ffProjCode("04-brca_DNAm_figures.R"))
#############################################################################

# source(ffProjCode("19-permTestBRCA_ATAC_PCA.R"))
# source(ffProjCode("06_2-brca_ATAC_vis")) # self-contained (doesn't run main visualization pipeline)

#############################################################################

# source(ffProjCode("19-permTestMOFACLLDNAm.R"))
# source(ffProjCode("09-MOFA_visualization.R")) # self-contained (loads data)
# source(ffProjCode("09-MOFA_figures.R"))

#############################################################################

source(ffProjCode("22-otherCancerAttrCorrelation.R")) # KIRC analysis
# KIRC, haven't run comprehensive visualization
source(ffProjCode("23-survivalAnalysis")) # KIRC, self-cont?
source(ffProjCode("22-KIRC_figures.R"))

#############################################################################
