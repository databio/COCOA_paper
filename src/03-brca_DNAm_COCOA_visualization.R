# 
# PCsToAnnotate_cPCH = c("PC1m4", "PC1p3", PCSTOANNOTATE)

source(paste0(Sys.getenv("CODE"),"COCOA_paper/src/00-init.R"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
########################################
source(paste0(Sys.getenv("CODE"), "/COCOA_paper/src/load_process_regions_brca.R"))

#################################################################

# parameters for COCOA_vis_pipeline.R
plotSubdir = "03_brca_COCOA_Vis"
# a cache is created in the script
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/")) 
simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "bigSharedC", reload = TRUE)
# screen out patients without ER/PGR status
methylData = bigSharedC[["methylProp"]][, 
                                          colnames(bigSharedC[["methylProp"]]) %in% patientMetadata[, subject_ID]] 
inputID = "sharedC"

# GRList # from load_process_regions pipeline
coordinateDT = bigSharedC$coordinates
allMPCAString = "allMPCA_657"
simpleCache(allMPCAString, assignToVariable = "mPCA", reload = TRUE)
loadingMat = mPCA$rotation
# use rsEnString to specify?
simpleCache("rsEnrichment_657", assignToVariable = "rsEnrichment", reload = TRUE)
simpleCache("rsEnrichmentTop10_657", assignToVariable = "rsEnrichmentTop10", reload = TRUE)
# TODO make sure GRList and rsEnrichment are both in the same order/with same data
names(GRList) <- paste0(rsEnrichment$rsName, " : ", rsEnrichment$rsDescription)
GRList = GRList[!is.na(rsEnrichment$PC1)]
rsEnrichment=rsEnrichment[!is.na(rsEnrichment$PC1), ]
rsEnSortedInd= rsRankingIndex(rsEnrichment = rsEnrichment, PCsToAnnotate = paste0("PC", 1:10))

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

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/COCOA_vis_pipeline.R"))

###############################################################################

# plotting correlation between a PC and the "PC-subset" score 
# derived from only loading values of CpGs within a certain region set 
# do this for top few region sets for PC1
PCofInterest = paste0("PC", 1:5)
for (j in seq_along(PCofInterest)) {
grDevices::pdf(paste0(Sys.getenv("PLOTS"), 
                      plotSubdir, "/subsetPC_PC_Scatter_", 
                      PCofInterest[j], addUnderscore(inputID, side="left"), ".pdf"))
for (i in 1:10) {
    # returns a "PC-subset" score for each sample
    subScores = pcFromSubset(regionSet = GRList[as.numeric(rsEnSortedInd[i, PCofInterest[j], with=FALSE])], 
                    pca = mPCA, 
                    methylData = methylData, 
                    coordinateDT = coordinateDT, 
                    PCofInterest = PCofInterest[j],
                    returnCor = FALSE)    
    
    plot(x = subScores, y = mPCA$x[, PCofInterest[j]], xlab= "Scores from only CpGs in this region set",
         ylab = PCofInterest[j], main = paste0(rsEnrichment$rsName[as.numeric(rsEnSortedInd[i, PCofInterest[j], with=FALSE])], 
                                            " : ", rsEnrichment$rsDescription[as.numeric(rsEnSortedInd[i, PCofInterest[j], with=FALSE])]))   

}
dev.off()

}


############################## just generating plots of PC vs. PC score from subset of CpGs 
rawScoreRSEn = get(load(paste0(getOption("RCACHE.DIR"), "/old_RCache/raw_loading_average/rsEnrichment_657.RData")))

if (!dir.exists(paste0(Sys.getenv("PLOTS"), plotSubdir, "/raw/"))) {
    dir.create(paste0(Sys.getenv("PLOTS"), plotSubdir, "/raw/"), recursive = TRUE)
}

rawScoreInd = rsRankingIndex(rawScoreRSEn, PCsToAnnotate = "PC1")
PCofInterest = paste0("PC", 1)
for (j in seq_along(PCofInterest)) {
    grDevices::pdf(paste0(Sys.getenv("PLOTS"), 
                          plotSubdir, "/raw/subsetPC_PC_Scatter_", 
                          PCofInterest[j], addUnderscore(inputID, side="left"), ".pdf"))
    for (i in 1:10) {
        # returns a "PC-subset" score for each sample
        subScores = pcFromSubset(regionSet = GRList[as.numeric(rawScoreInd[i, PCofInterest[j], with=FALSE])], 
                                 pca = mPCA, 
                                 methylData = methylData, 
                                 coordinateDT = coordinateDT, 
                                 PCofInterest = PCofInterest[j],
                                 returnCor = FALSE)    
        
        plot(x = subScores, y = mPCA$x[, PCofInterest[j]], xlab= "Scores from only CpGs in this region set",
             ylab = PCofInterest[j], main = paste0(rsEnrichment$rsName[as.numeric(rawScoreInd[i, PCofInterest[j], with=FALSE])], 
                                                   " : ", rsEnrichment$rsDescription[as.numeric(rawScoreInd[i, PCofInterest[j], with=FALSE])]))   
        
    }
    dev.off()
    
}


####### generating meta-region plots with only top 4 region sets