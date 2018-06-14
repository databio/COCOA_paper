# Plots for ewing_PCRSA_visualization

source(paste0(Sys.getenv("CODE"),"pcrsa_method_paper/src/00-init.R"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "ews_patients/analysis/plots/"))
########################################
source(paste0(Sys.getenv("CODE"), "/pcrsa_method_paper/src/load_process_regions_brca.R"))

#################################################################

# parameters for PCRSA_vis_pipeline.R
plotSubdir = "53_Ewing_PCRSA_Vis"
# a cache is created in the script
setCacheDir(paste0(Sys.getenv("PROCESSED"), "ews_patients/RCache/")) 
simpleCache("bigSharedC_pQC", assignToVariable = "bigSharedC", reload = TRUE)

# GRList # from load_process_regions pipeline
coordinateDT = bigSharedC$coordinates
methylData = bigSharedC$methylProp
allMPCAString = "allMPCA2"
simpleCache(allMPCAString, assignToVariable = "mPCA", reload = TRUE)
loadingMat = mPCA$rotation
# use rsEnString to specify?
simpleCache("rsEnrichment") # , assignToVariable = "rsEnrichment")
simpleCache("rsEnrichmentTop10") # , assignToVariable = "rsEnrichmentTop10")
# TODO make sure GRList and rsEnrichment are both in the same order/with same data
names(GRList) <- paste0(rsEnrichment$rsName, " : ", rsEnrichment$rsDescription)
GRList = GRList[!is.na(rsEnrichment$PC1)]
rsEnrichment=rsEnrichment[!is.na(rsEnrichment$PC1), ]
rsEnSortedInd= rsRankingIndex(rsEnrichment = rsEnrichment, PCsToAnnotate = paste0("PC", 1:10))
## "region set Overlapping Cytosine Proportion" (rsOLCP)
## proportion of cytosines from region set that are shared with other region set
rsEnrichmentTop10 = rsEnrichmentTop10[!is.na(rsEnrichmentTop10$PC1), ]
PCSTOANNOTATE = paste0("PC", 1:10)


### plots that will be created and script specific parameters for them 
# "comparePCHeatmap"
PCsToAnnotate_cPCH = PCSTOANNOTATE
# "methylAlongPC"
topRSToPlotNum = 15
PCsToAnnotate_mAPC = PCSTOANNOTATE[1:5]
# "regionQuantileByPC"
PCsToAnnotate_rQBPC = PCSTOANNOTATE
topRSInd_rQBPC = unique(unlist(rsEnSortedInd[1:15, ])) # get top region sets from each PC
# pcFromSubset Correlation Heatmap
PCsToAnnotate_pcFSCH = PCSTOANNOTATE
topRSInd_pcFSCH = unique(unlist(rsEnSortedInd[1:15, ])) # get top region sets from each PC
## "region set Overlapping Cytosine Proportion" (rsOLCP)
## proportion of cytosines from region set that are shared with other region set
topRSInd_rsOLCP = unique(unlist(rsEnSortedInd[1:15, ]))

# the pipeline
source("PCRSA_vis_pipeline.R") 