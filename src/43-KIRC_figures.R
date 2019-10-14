
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(cowplot)
library(ROCR)

# loadProcessKIRCMethyl()


################################################################################
# figures for KIRC supervised example

# simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), assignToVariable = "rsPermScores")
simpleCache(paste0("rsScores_",  dataID, "_", variationMetric), assignToVariable = "realRSScores")
p = plotRSConcentration(rsScores = realRSScores, scoreColName = "cancerStage", pattern = "EZH2|SUZ12")
ggsave(filename = ffPlot(paste0(plotSubdir, "cancerStage", 
                                "EZH2_SUZ12_", dataID, "_", variationMetric, ".svg")), plot = p, device = "svg")


########
# new figure that has bar for each region set and shows actual scores
#https://bioinfo.iric.ca/introduction-to-cowplot/
plot(1:nrow(realRSScores), realRSScores[order(realRSScores$cancerStage, decreasing = TRUE), ]$cancerStage)



a = plotAnnoScoreDist2(rsScores = realRSScores, colsToPlot = "cancerStage", pattern = c("EZH2|SUZ12", "H3K27me3"))
a + ggtitle("Association between DNA methylation and cancer stage ")
plotAnnoScoreDist(rsScores = realRSScores, colsToPlot = "cancerStage", pattern = c("EZH2|SUZ12"))
plotAnnoScoreDist(rsScores = realRSScores, colsToPlot = "cancerStage", pattern = c("nonefound"))



# either color dots for EZH2/Suz12 or have color bar beside plot
# make region sets that I will be pulling out to plot separately a different color also (e.g. JunD)?
# have a color bar beside plot marking whether each RS was statistically significant or not?

########
simpleCache(paste0("rsPermScores_",  nPerm, "Perm_", variationMetric, "_", dataID), 
            assignToVariable = "rsPermScores")

# include null distributions
permResultsMat = do.call(cbind, lapply(rsPermScores, function(x) x$cancerStage))

# make a ROC curve plot for EZH2/Suz12
# pred = realRSScores$cancerStage / max(realRSScores$cancerStage, na.rm = TRUE)
pred = cbind(realRSScores$cancerStage, permResultsMat)
rocPreds = ROCR::prediction(pred, labels = matrix(grepl(pattern = "EZH2|SUZ12", 
                                        x = realRSScores$rsName, 
                                        ignore.case = TRUE), nrow=nrow(pred), ncol=ncol(pred)))
testAUC = ROCR::performance(rocPreds, measure="auc")@y.values[[1]]
allTestAUC = unlist(ROCR::performance(rocPreds, measure="auc")@y.values)
testAUC
perf = ROCR::performance(rocPreds, measure = "tpr", x.measure = "fpr")
plot(perf)
title(main= "ROC curve for EZH2 and SUZ12")

