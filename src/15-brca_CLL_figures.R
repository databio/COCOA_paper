# library(projectInit)

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))

set.seed(1234)
plotSubdir = "15_brca_CLL_figures/"
scriptID = "15-brcaCLLFigures"
createPlotSubdir(plotSubdir = plotSubdir)

###########################################################
# reading in the region sets
# load LOLA database

loadGRList(genomeV = "hg38")
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


##############################################################################
# Fig. 4, MOFA CLL
# http://msb.embopress.org/content/14/6/e8124
# factors 1 (strong), 7 and 9 were associated with DNA methylation
# factor 1: cell type/differentiation, factor 7: chemo-immunotherapy treatment prior to sample collection
# factor 7: del17p, TP53 mutations, methylation of oncogenes

loadMOFAData(methylMat = TRUE, signalCoord=TRUE, latentFactors = TRUE,
             cllMultiOmics=TRUE)
simpleCache("mofaPermZScoresCor", assignToVariable = "rsZScores")
simpleCache("rsScore_Cor_CLL196MOFA", assignToVariable = "rsScores")

mutDF = cllMultiOmics$Mutation
# same order
latentFactors = latentFactors[colnames(mutDF), ]
IGHV = mutDF["IGHV", ]
IGHV[IGHV == 1] = "mutated"
IGHV[IGHV == 0] = "unmutated"
IGHV = as.factor(IGHV)
latentFactors = data.frame(latentFactors, IGHV, KLHL6=mutDF["KLHL6", ], BRAF=factor(mutDF["BRAF", ])) 
                           # CD56=cllMultiOmics$mRNA["ENSG00000149294", ],
                           # CD3=cllMultiOmics$mRNA["ENSG00000167286", ])
# CD19=cllMultiOmics$mRNA["ENSG00000177455", ]
# latentFactors$IGHV = factor(latentFactors$IGHV, levels = c("mutated", "unmutated"))

# make plot of latent factors (show that LF1 separates based on differentiation)
# it seems that the sign of latent factor 1 is switched either in the paper or the data package
# -1 * LF1 to make it the same as the paper
fPlot = ggplot(data = latentFactors, mapping = aes(x = -1*LF1, y = LF2)) + geom_point(aes(col=IGHV)) +
            xlab("Factor 1") + ylab("Factor 2") + theme(axis.title.x = element_text(size = 18),
                                                        axis.title.y = element_text(size = 18),
                                                        legend.text = element_text(size=18), 
                                                        legend.title = element_text(size=18))
fPlot
ggsave(filename = ffPlot(paste0(plotSubdir, "mofaFactorLF1LF2.svg")), 
       plot = fPlot, device = "svg")
cor.test(x = latentFactors$CD3, y = latentFactors$LF9)



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

#################### plot raw data in top regions for LF1
simpleCache("inferredMethylWeightsMOFA", assignToVariable = "featureLFCor")
# make sure they are in the same order/have same CpGs
featureLFCor = featureLFCor[row.names(methylMat), ]
loadGRList(genomeV = "hg19")
dataID = "cll196"
lfCols= paste0("LF", c(1:3, 5:7, 9))
cpgCov = rsScores$cytosine_coverage
topRSInd = rsRankingIndex(rsScores = rsScores[cpgCov >= 200, ], signalCol = lfCols)
View(rsScores[order(rsScores$LF1, decreasing = TRUE), ])
topLF1Ind = topRSInd$LF1[1:20]



# top few region sets for LF1
# arbitrarily selecting region sets that covered at least 200 CpGs
topLF1RS = GRList[cpgCov >= 200][topLF1Ind]
topLF1RSNames = rsName[cpgCov >= 200][topLF1Ind]
topLF1RSDes = rsDescription[cpgCov >= 200][topLF1Ind]
    
for (i in seq_along(lfCols)) {
    theseTopRSInd = as.numeric(topRSInd[1:20, lfCols[i]])
    theseTopRS = GRList[cpgCov >= 200][theseTopRSInd]
    theseTopRSNames = rsName[cpgCov >= 200][theseTopRSInd]
    theseTopRSDes = rsDescription[cpgCov >= 200][theseTopRSInd]
    
    # # top region sets for this PC
    # rsInd = as.numeric(as.matrix(rsEnSortedInd[1:topRSToPlotNum, lfCols[i]])) # original index
    
    grDevices::pdf(ffPlot(paste0(plotSubdir, "regionMethylHeatmaps", lfCols[i], dataID, ".pdf")), width = 11, height = 8.5)
    
    # heatmap
    for (j in seq_along(theseTopRS)) {
        print(signalAlongPC(genomicSignal=methylMat,
                      signalCoord=signalCoord,
                      sampleScores=latentFactorMat,
                      regionSet=theseTopRS[[j]], orderByCol=lfCols[i],
                      topXVariables=50,
                      variableScores = abs(as.numeric(featureLFCor[, lfCols[i]])),
                      cluster_columns = TRUE, column_title = theseTopRSNames[j]))
    }
    
    # draw(Heatmap(matrix = methylData[1:1000, 1:10]))
    # plot(methylData[1:1000, 1])
    
    dev.off()
}

signalAlongPC(genomicSignal=methylMat,
              signalCoord=signalCoord,
              sampleScores=latentFactorMat,
              regionSet=topLF1RS[["Gm12878_WE.bed"]], orderByCol=lfCols[i],
              topXVariables=50,
              variableScores = abs(as.numeric(featureLFCor[, lfCols[i]])),
              cluster_columns = TRUE, column_title = "Individual cytosines")


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




