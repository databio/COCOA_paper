# library(projectInit)

# project.init(codeRoot = paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/"), dataDir = paste0(Sys.getenv("PROCESSED"), "brca_PCA/"))
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
#patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)
plotSubdir = "15_brca_CLL_figures/"

###########################################################
# reading in the region sets
# load LOLA database
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))
names(GRList) = rsName

#################################################################

dataID = "657" # 657 patients with both ER and PGR info in metadata, 692 total
allMPCAString = "allMPCA_657" #  "allMPCA_657"
top10MPCAString = "top10MPCA_657"

simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaSharedC", reload = TRUE)


simpleCache("rsEnrichmentRSMean_657", assignToVariable = "rsScores")

# get top region sets for DNA methylation results
topScoreResultsInd = rsRankingIndex(rsScores, "PC1")$PC1[1:20]
topGRList = GRList[rsScores$rsName[topScoreResultsInd]]

# loading PCA and combining components that could separate ER=/-
# for rsEnrichment, PCs 1 and 4 could separate ER+/-
simpleCache(allMPCAString, assignToVariable = "mPCA")
brcaLoadings = mPCA$rotation
brcaCoord = brcaSharedC[["coordinates"]]

#############################################################################
# Fig. 2a, PCA plot, DNA methylation

colorByCols = "ER_status"
pcaWithAnno = cbind(mPCA$x, patientMetadata[row.names(mPCA$x) ,])
multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                       plotCols = c("PC1", "PC4"), 
                                       colorByCols=colorByCols)
ggplot2::ggsave(filename=paste0(Sys.getenv("PLOTS"), "/allMPCA_PCA_Plots/multiColorPCAPlots_allMPCA_", 
                                PCsToPlot[i], "x", PCsToPlot[j], 
                                ".pdf"), plot = multiColorPCAPlots, device = "pdf",
                limitsize=FALSE)


##################
# Fig 2. BRCA DNA methylation 
simpleCache("rsEnrichmentRSMean_657", assignToVariable = "rsScores")

# summary figure of COCOA BRCA results: ER set relative ranking among region sets
esrConcentrationPlot = plotRSConcentration(rsScores, scoreColName=paste0("PC", 1), 
                            colsToSearch = c("rsName", "rsDescription"), 
                            pattern= "esr|eralpha|eraa") + ggtitle("Estrogen receptor region sets") + 
    theme(axis.title.x = element_text(size = 20), 
                             axis.title.y = element_text(size=20), 
          axis.text.x = element_text(size=20), 
          axis.text.y = element_text(size=20), title = element_text(size = 20)) + scale_y_continuous(breaks = seq(from=0, to=10, by=2))
ggsave(filename = ffPlot(paste0(plotSubdir, "erRegionSetsDNAm.svg")), plot = esrConcentrationPlot, device = "svg")

erRelated = plotRSConcentration(rsScores, scoreColName=paste0("PC", 1), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|eralpha|gata3|foxa1|h3r17") + ggtitle("Estrogen receptor-related region sets")
ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "ER_related_PC1.svg"), plot = erRelated, device = "svg")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|eralpha|gata3|foxa1|h3r17")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "ezh2|suz12")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "pol2")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k9")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k27me")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k9|h3k27me")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k4me3")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k36")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "k562")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "mcf7|mcf-7")
###################
# Fig. 2c, meta region loading profile plots, DNA methylation BRCA

PCsToAnnotate = paste0("PC", 1:4)
wideGRList <- lapply(topGRList, resize, width=14000, fix="center")
loadProfile <- lapply(wideGRList, function(x) getLoadingProfile(loadingMat=brcaLoadings,
                                                                signalCoord=brcaCoord,
                                                                regionSet=x, 
                                                                PCsToAnnotate=PCsToAnnotate,
                                                                binNum=21))

# average loading value from each PC to normalize so PCs can be compared with each other
if (is.numeric(brcaLoadings[, PCsToAnnotate])) {
    avLoad = mean(abs(brcaLoadings[, PCsToAnnotate]))
} else {
    avLoad <- apply(X=brcaLoadings[, PCsToAnnotate], 
                    MARGIN=2, 
                    FUN=function(x) mean(abs(x)))
}

# normalize
loadProfile <- lapply(loadProfile, 
                      FUN=function(x) as.data.frame(mapply(FUN = function(y, z) x[, y] - z, 
                                                           y=PCsToAnnotate, z=avLoad)))
binID = 1:nrow(loadProfile[[1]])
loadProfile <- lapply(loadProfile, FUN=function(x) cbind(binID, x))
# for the plot scale
maxVal <- max(sapply(loadProfile, FUN=function(x) max(x[, PCsToAnnotate])))
minVal <- min(sapply(loadProfile, FUN=function(x) min(x[, PCsToAnnotate])))
# convert to long format for plots
loadProfile <- lapply(X=loadProfile, FUN=function(x) tidyr::gather(data=x, key="PC", value="loading_value", PCsToAnnotate))
loadProfile <- lapply(loadProfile, 
                      function(x){x$PC <- factor(x$PC, levels=PCsToAnnotate); return(x)})
wrapper <- function(x, ...) paste(strwrap(x, ...), collapse="\n") 

# manually get region set name
regionSetNames = c("Gata3", "H3R17me2", "ESR1", "Gata3", 
                   "FoxA1", "ESR1", "ESR1", "ESR1", 
                   "FoxA1", "ESR1", "AR", "AR", 
                   "Znf217", "TCF7L2", "ESR1", "ESR1", 
                   "JunD", "ESR1", "FoxA1", "H3R17me2")

profilePList <- list()
for (i in seq_along(loadProfile)) {
    
    thisRS <- loadProfile[[i]]
    
    profilePList[[i]] <- ggplot(data=thisRS, 
                                mapping=aes(x=binID , y=loading_value)) + 
        geom_line() + ylim(c(minVal, maxVal)) + facet_wrap(facets="PC") + 
        ggtitle(label=wrapper(regionSetNames[i], width=30)) + 
        xlab("Genome around region set, 14 kb") + 
        ylab("Normalized loading value") + 
        theme(panel.grid.major.x=element_blank(), 
              panel.grid.minor.x=element_blank(), 
              axis.text.x=element_blank(), 
              axis.ticks.x=element_blank())
    ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "metaRegionPlots/", regionSetNames[i], "_", i), plot = profilePList[[i]], device = "pdf")
}
profilePList[[1]]
ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, regionSetName) )

profilePList[[2]]
profilePList[[3]]
profilePList[[4]]


########
### figures for PC4, Supplementary?
View(rsScores[order(rsScores$PC4, decreasing = TRUE), ])
#
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "mcf7|mcf-7")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k9|h3k27me|suz12|ezh2")
# rsConcentration is only plotting first 25 in bin 1 (-25 to +25 instead of 0-50)?

# meta-region loading profile

###############################################################################
# Fig. 3, BRCA ATAC

# PCA plots


############################
# Fig. 4, MOFA CLL
# http://msb.embopress.org/content/14/6/e8124
# factors 1 (strong), 7 and 9 were associated with DNA methylation
# factor 1: cell type/differentiation, factor 7: chemo-immunotherapy treatment prior to sample collection
# factor 7: del17p, TP53 mutations, methylation of oncogenes
simpleCache("rsScore_Cor_CLL196MOFA", assignToVariable = "rsScores")
# View(rsScores[order(rsScores$LF1, decreasing=TRUE), ])
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                                colsToSearch = c("rsName", "rsDescription"), 
                                pattern= "K562")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|eralpha|gata3|foxa1|h3r17")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "GM12878|GM18526|GM12891|GM10847|K562|leukemia|leukaemia|lymphoma")
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

############################################################################
# region set results distribution plot (base R plotting)



# plotRSConcentration <- function(rsScores, scoreColName="PC1", 
#                                 colsToSearch = c("rsName", "rsDescription"), 
#                                 pattern, breaks) {
#     
#     rsScores = as.data.frame(rsScores)
#     rsScores = rsScores[order(rsScores[, scoreColName], decreasing = TRUE), ]
#     
#     rsInd = rep(FALSE, nrow(rsScores))
#     for (i in seq_along(colsToSearch)) {
#         rsInd = rsInd | grepl(pattern = pattern, x = rsScores[, colsToSearch[i]], ignore.case = TRUE)
#     }
#     
#     rsInd = which(rsInd)
#     hist(rsInd, breaks = breaks)
#     
# }
# 
# plotRSConcentration(rsScores=rsScores, scoreColName = "PC5", pattern = "esr1|eralpha", breaks=seq(0, 2000, by=200))
# plotRSConcentration(rsScores=rsScores, 
#                     scoreColName = "PC1", 
#                     pattern = "stat", 
#                     breaks=seq(0, 2400, by=200))




