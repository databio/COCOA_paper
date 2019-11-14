# 

source(paste0(Sys.getenv("CODE"),"COCOA_paper/src/00-init.R"))
#####################################################################

# parameters for COCOA_vis_pipeline.R
plotSubdir = "09_MOFA_Vis/"
# dataID = "CLL196MOFA"
# variationMetric = "cov"
# nPerm = 300
rsScoreCacheName = paste0("rsScore_", dataID, "_", variationMetric)
inputID = paste0(nPerm, "Perm_", dataID, "_", variationMetric) # used by COCOA vis script

#################################################################
# load hg19 database
loadGRList(genomeV = "hg19")

#################################################################
# load MOFA CLL data
cllMethyl = prepareCLLMethyl()
methylData = cllMethyl$methylProp
methCoord = cllMethyl$methylCoord

# GRList # from load_process_regions pipeline
coordinateDT = cllMethyl$methylCoord
coordinateDT = COCOA:::grToDt(coordinateDT)

# simpleCache(allMPCAString, assignToVariable = "mPCA", reload = TRUE)
simpleCache(paste0("inferredMethylWeightsMOFA", "_", variationMetric), assignToVariable = "loadingMat")
simpleCache("cllMOFAFactors", assignToVariable = "latentFactors")
# use rsEnString to specify?
simpleCache(rsScoreCacheName, assignToVariable = "rsEnrichment", reload = TRUE)
rsEnrichment = cbind(rsEnrichment, rsCollection)

# screen out region sets with low coverage of the data
lowCov = rsEnrichment$signalCoverage < 100
rsEnrichment = rsEnrichment[!lowCov, ]
GRList = GRList[!lowCov]

# screen out roadmap epigenomics sets since they aren't very interpretable
keepInd = screenOutRoadmap(rsScores = rsEnrichment, 
                 rsCollection = rsEnrichment$rsCollection, patternSearchCol="rsName", keepPattern="H3K4me1") 
rsEnrichment = rsEnrichment[keepInd, ]
GRList = GRList[keepInd]

# the latent factors
mPCA = list()
mPCA$x = latentFactors
mPCA$center = rep(0, nrow(methylData))
mPCA$rotation = loadingMat
if (!all(row.names(latentFactors) == colnames(methylData))) {
    stop("samples are not ordered consistently")
}

# TODO make sure GRList and rsEnrichment are both in the same order/with same data
names(GRList) <- paste0(rsEnrichment$rsName, " : ", rsEnrichment$rsDescription)
GRList = GRList[!is.na(rsEnrichment$LF1)]
rsEnrichment=rsEnrichment[!is.na(rsEnrichment$LF1), ]
rsEnSortedInd= rsRankingIndex(rsScores = rsEnrichment, signalCol = paste0("LF", 1:10))

################################################################################
# PCSTOANNOTATE = paste0("LF", c(1:3, 5:7, 9))
PCSTOANNOTATE = paste0("LF", 1:10)

### plots that will be created and script specific parameters for them 
# "comparePCHeatmap"
PCsToAnnotate_cPCH = PCSTOANNOTATE
# "methylAlongPC"
topRSToPlotNum = 15
PCsToAnnotate_mAPC = PCSTOANNOTATE
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

###################

# the pipeline
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/COCOA_vis_pipeline.R")) 

################################################################################
# # creating figures for presentation that only include a few PCs
# plotSubdir = 
# inputID =  
# PCSTOANNOTATE = paste0("PC", 1:4)
# 
# ### plots that will be created and script specific parameters for them 
# # "comparePCHeatmap"
# PCsToAnnotate_cPCH = PCSTOANNOTATE
# # "methylAlongPC"
# topRSToPlotNum = 15
# PCsToAnnotate_mAPC = PCSTOANNOTATE[1:5]
# # "regionQuantileByPC"
# PCsToAnnotate_rQBPC = PCSTOANNOTATE
# topRSInd_rQBPC = unique(unlist(rsEnSortedInd[1:15, ])) # get top region sets from each PC
# # pcFromSubset Correlation Heatmap
# PCsToAnnotate_pcFSCH = PCSTOANNOTATE
# topRSInd_pcFSCH = unique(unlist(rsEnSortedInd[1:10, c("PC1", "PC4"), with=FALSE])) # get top region sets from each 
# ## "region set Overlapping Cytosine Proportion" (rsOLCP)
# ## proportion of cytosines from region set that are shared with other region set
# topRSInd_rsOLCP = unique(unlist(rsEnSortedInd[1:15, c("PC1", "PC2", "PC3", "PC4"), with=FALSE]))
# ## "meta region loading profiles" (mrLP)
# topRSInd_mrLP = unique(unlist(rsEnSortedInd[1:15, c(paste0("PC", 1:4)), with=FALSE]))
# PCsToAnnotate_mrLP = PCSTOANNOTATE
# 
# ##################
# source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/COCOA_vis_pipeline.R"))
# 
# ###############################################################################
# old exploratory visualization

# View(rsScores[order(rsScores$LF1, decreasing=TRUE), ])
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "K562")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "GM12878|GM18526|GM12891|GM10847|K562|leukemia|leukaemia|lymphoma")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|eralpha|gata3|foxa1|h3r17")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k9")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k4me1")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k4me3")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k36")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k27me")
# transcription factors
# associated with immune latent factors
plotRSConcentration(rsScores[rsScores$region_coverage >= 100, ], scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "nfkb")
# 
plotRSConcentration(rsScores[rsScores$region_coverage >= 100, ], scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr1")

#############################################################################