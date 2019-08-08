
source(paste0(Sys.getenv("CODE"),"COCOA_paper/src/00-init.R"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))

# parameters for COCOA_vis_pipeline.R
plotSubdir = "03_brca_COCOA_Vis"
inputID = "sharedC"

# a cache is created in the script
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/")) 
########################################
source(paste0(Sys.getenv("CODE"), "/COCOA_paper/src/load_process_regions_brca.R"))

#################################################################


simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "bigSharedC", reload = TRUE)
# screen out patients without ER/PGR status
methylData = bigSharedC[["methylProp"]][, 
                                        colnames(bigSharedC[["methylProp"]]) %in% patientMetadata[, subject_ID]] 


# GRList # from load_process_regions pipeline
coordinateDT = bigSharedC$coordinates
allMPCAString = "allMPCA_657"
simpleCache(allMPCAString, assignToVariable = "mPCA", reload = TRUE)
loadingMat = mPCA$rotation
# use rsEnString to specify?
simpleCache("rsEnrichmentRSMean_657", assignToVariable = "rsEnrichment", reload = TRUE)
simpleCache("rsEnrichmentTop10_657", assignToVariable = "rsEnrichmentTop10", reload = TRUE)
# TODO make sure GRList and rsEnrichment are both in the same order/with same data
names(GRList) <- paste0(rsEnrichment$rsName, " : ", rsEnrichment$rsDescription)
GRList = GRList[!is.na(rsEnrichment$PC1)]
rsEnrichment=rsEnrichment[!is.na(rsEnrichment$PC1), ]
rsEnSortedInd= rsRankingIndex(rsScores = rsEnrichment, PCsToAnnotate = paste0("PC", 1:10))

rsEnrichmentTop10 = rsEnrichmentTop10[!is.na(rsEnrichmentTop10$PC1), ]
PCSTOANNOTATE = paste0("PC", 1:10)

#################################################################################
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
topRSInd_pcFSCH = unique(unlist(rsEnSortedInd[1:15, ])) # get top region sets from each 
## "region set Overlapping Cytosine Proportion" (rsOLCP)
## proportion of cytosines from region set that are shared with other region set
topRSInd_rsOLCP = unique(unlist(rsEnSortedInd[1:10, ]))
## "meta region loading profiles" (mrLP)
topRSInd_mrLP = unique(unlist(rsEnSortedInd[1:10, ]))
PCsToAnnotate_mrLP = PCSTOANNOTATE

# the pipeline
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/COCOA_vis_pipeline.R")) 


# creating figures for presentation that only include a few PCs
plotSubdir = "03_brca_pres_figures"
inputID = "sharedC_figures"

PCSTOANNOTATE = paste0("PC", 1:4)

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
topRSInd_pcFSCH = unique(unlist(rsEnSortedInd[1:10, c("PC1", "PC4"), with=FALSE])) # get top region sets from each 
## "region set Overlapping Cytosine Proportion" (rsOLCP)
## proportion of cytosines from region set that are shared with other region set
topRSInd_rsOLCP = unique(unlist(rsEnSortedInd[1:15, c("PC1", "PC2", "PC3", "PC4"), with=FALSE]))
## "meta region loading profiles" (mrLP)
topRSInd_mrLP = unique(unlist(rsEnSortedInd[1:15, c(paste0("PC", 1:4)), with=FALSE]))
PCsToAnnotate_mrLP = PCSTOANNOTATE

#################################################################################
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/COCOA_vis_pipeline.R"))

###############################################################################