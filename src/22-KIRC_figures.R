
################################################################################
# figures for KIRC supervised example
dataID = "kircMethyl"
simpleCache(paste0("rsScores_", dataID, "Cor"), assignToVariable = "realRSScores")
p = plotRSConcentration(rsScores = realRSScores, scoreColName = "cancerStage", pattern = "EZH2|SUZ12")
ggsave(filename = ffPlot(paste0(plotSubdir, "cancerStage", "EZH2_SUZ12_", dataID, ".svg")), plot = p, device = "svg")
