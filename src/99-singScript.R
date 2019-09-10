# script to run on singularity on cluster

.libPaths(c("/scratch/jtl2hk/3.6", .libPaths()))

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

#############################################################################
# source(ffProjCode("19-permTestBRCADNAmPCA.R"))
# source(ffProjCode("03-brca_DNAm_COCOA_visualization")) # self-contained
# source(ffProjCode("04-brca_DNAm_figures.R"))
#############################################################################

source(ffProjCode("19-permTestBRCA_ATAC_PCA.R"))
source(ffProjCode("06_2-brca_ATAC_vis")) # self-contained (doesn't run main visualization pipeline)

#############################################################################




#############################################################################



#############################################################################
