


########################################

source(paste0(Sys.getenv("CODE"), "/pcrsa_method_paper/src/load_process_regions_brca.R"))

#################################################################

# parameters for PCRSA_vis_pipeline.R
plotSubdir = "53_Ewing_PCRSA_Vis"
setCacheDir(paste0(Sys.getenv("PROCESSED"), "ews_patients/RCache/")) 
rsEnSortedInd
# GRList # from load_process_regions pipeline 
coordinateDT = brcaMList$coordinates
loadingMat = allMPCA$rotation
methylData = brcaMList$methylProp
allMPCAString 
simpleCache(allMPCAString, assignToVariable = "mPCA")
# use rsEnString to specify?
simpleCache("rsEnrichment_657", assignToVariable = "rsEnrichment")
simpleCache("rsEnrichmentTop10_657", assignToVariable = "rsEnrichmentTop10")

PCSTOANNOTATE = paste0("PC", 1:10)




### plots that will be created and script specific parameters for them 
# "comparePCHeatmap"
PCsToAnnotate_cPCH = PCSTOANNOTATE
# "methylAlongPC"
PCsToAnnotate_mAPC = PCSTOANNOTATE[1:5]
# "regionQuantileByPC"
PCsToAnnotate_rQBPC = c("PC1m4", "PC1p3", PCSTOANNOTATE)
topRSInd_rQBPC = unique(unlist(rsEnSortedInd[1:10, ])) # get top 10 region sets from each PC
# pcFromSubset Correlation Heatmap
PCsToAnnotate_pcFSCHM = PCSTOANNOTATE
topRSInd_pcFSCHM = unique(unlist(rsEnSortedInd[1:10, ])) # get top 10 region sets from each PC

# the pipeline
source("PCRSA_vis_pipeline.R") 