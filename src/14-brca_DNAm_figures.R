source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))

set.seed(1234)
plotSubdir = "99_paper_figures/"
scriptID = "99-paperFigures"
createPlotSubdir(plotSubdir = plotSubdir)
.analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)

###########################################################
# reading in the region sets
# load LOLA database

loadGRList(genomeV = "hg38")

#################################################################

allMPCAString = "allMPCA_657" #  "allMPCA_657"

simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaSharedC", reload = TRUE)

simpleCache(paste0("rsScores_", dataID, "_", variationMetric), assignToVariable = "rsScores")

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
ggplot2::ggsave(filename=ffPlot(paste0(plotSubdir,"/allMPCA_PCA_Plots/multiColorPCAPlots_allMPCA_", 
                                PCsToPlot[i], "x", PCsToPlot[j], 
                                ".pdf")), plot = multiColorPCAPlots, device = "pdf",
                limitsize=FALSE)

##################
# Fig 2. BRCA DNA methylation 
# simpleCache(paste0("rsScores_", dataID, "_", variationMetric), assignToVariable = "rsScores")

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
####################
# annoScoreDist
a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = "PC1", pattern = "esr1|eralpha")
a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = "PC1", pattern = c("gata3|foxa1|H3R17me2", "esr1|eralpha"))
a
a
###################
# Fig. 2c, meta region loading profile plots, DNA methylation BRCA

PCsToAnnotate = paste0("PC", 1:4)
# # manually get region set name
# regionSetNames = c("Gata3", "H3R17me2", "ESR1", "Gata3", 
#                    "FoxA1", "ESR1", "ESR1", "ESR1", 
#                    "FoxA1", "ESR1", "AR", "AR", 
#                    "Znf217", "TCF7L2", "ESR1", "ESR1", 
#                    "JunD", "ESR1", "FoxA1", "H3R17me2")
regionSetNames = names(topGRList)
wideGRList <- lapply(topGRList, resize, width=14000, fix="center")
mrProfileList  = makeMetaRegionPlots(signal=brcaLoadings, signalCoord=brcaCoord,
                    GRList=wideGRList, rsNames=regionSetNames, 
                    signalCol=PCsToAnnotate, binNum=21, aggrMethod="default")

# for (i in seq_along(loadProfile)) {
#     ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "metaRegionPlots/", regionSetNames[i], "_", i), plot = profilePList[[i]], device = "pdf")
# }
ggsave(plot = mrProfileList$grob, 
       filename = ffPlot(paste0(plotSubdir, "metaRegionPlots", .analysisID, ".pdf"), device = "pdf"))

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
# meta-region loading profile

########################################################################################
# ROC curve
simpleCache(paste0("rsPermScores_",  nPerm, "Perm_", variationMetric, "_", dataID), 
            assignToVariable = "rsPermScores")
realRSScores = rsScores
# include null distributions
permResultsMat = do.call(cbind, lapply(rsPermScores, function(x) x$PC1))

# make a ROC curve plot for EZH2/Suz12
# pred = realRSScores$cancerStage / max(realRSScores$cancerStage, na.rm = TRUE)
pred = cbind(realRSScores$PC1, permResultsMat)
rocPreds = ROCR::prediction(pred, labels = matrix(grepl(pattern = "esr1|eralpha", 
                                                        x = realRSScores$rsName, 
                                                        ignore.case = TRUE), nrow=nrow(pred), ncol=ncol(pred)))
testAUC = ROCR::performance(rocPreds, measure="auc")@y.values[[1]]
allTestAUC = unlist(ROCR::performance(rocPreds, measure="auc")@y.values)
testAUC
perf = ROCR::performance(rocPreds, measure = "tpr", x.measure = "fpr")
plot(perf)
title(main= "ROC curve for estrogen receptor")
