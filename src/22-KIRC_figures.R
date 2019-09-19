
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(cowplot)

################################################################################
# figures for KIRC supervised example
dataID = "kircMethyl"
simpleCache(paste0("rsScores_", dataID, "Cor"), assignToVariable = "realRSScores")
p = plotRSConcentration(rsScores = realRSScores, scoreColName = "cancerStage", pattern = "EZH2|SUZ12")
ggsave(filename = ffPlot(paste0(plotSubdir, "cancerStage", "EZH2_SUZ12_", dataID, ".svg")), plot = p, device = "svg")


########
# new figure that has bar for each region set and shows actual scores
#https://bioinfo.iric.ca/introduction-to-cowplot/
plot(1:nrow(realRSScores), realRSScores[order(realRSScores$cancerStage, decreasing = TRUE), ]$cancerStage)
realRSScores$rank = order(order(realRSScores$cancerStage, decreasing = TRUE))
realRSScores$polycomb = as.factor(as.numeric(grepl(pattern = "EZH2", 
                              x = realRSScores$rsName, 
                              ignore.case = TRUE) | grepl(pattern = "SUZ12", 
                                                          x = realRSScores$rsName, 
                                                          ignore.case = TRUE)))

rsCorP = ggplot(data = realRSScores, mapping = aes(x=rank, y=cancerStage)) + geom_point(aes(col=polycomb, alpha=polycomb)) 
rsCorP + coord_flip()
polycombStatusP = ggplot(data = realRSScores, mapping = aes(x=rank, y=as.numeric(polycomb)-1)) + geom_col(aes(col=polycomb)) 
polycombStatusP

gg_dist_g1 = polycombStatusP
gg_dist_g2 = rsCorP
gg_scatter = rsCorP

# gg_dist_g2 = gg_dist_g2 + coord_flip()

# Remove some duplicate axes
gg_dist_g1 = gg_dist_g1 + theme(axis.title.x=element_blank(),
                                axis.text=element_blank(),
                                axis.line=element_blank(),
                                axis.ticks=element_blank())

# Modify margin c(top, right, bottom, left) to reduce the distance between plots
#and align G1 density with the scatterplot
gg_dist_g1 = gg_dist_g1 + theme(plot.margin = unit(c(0.5, 0, 0, 0.5), "cm"))
gg_scatter = gg_scatter + theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm"))
#gg_dist_g2 = gg_dist_g2 + theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))

# Combine all plots together and crush graph density with rel_heights
first_col = plot_grid(gg_dist_g1, gg_scatter, ncol = 1, rel_heights = c(1, 20), align = "v")
#second_col = plot_grid(NULL, gg_dist_g2, ncol = 1, rel_heights = c(1, 3))
# perfect = plot_grid(first_col, second_col, ncol = 2, rel_widths = c(3, 1))
first_col

# either color dots for EZH2/Suz12 or have color bar beside plot
# make region sets that I will be pulling out to plot separately a different color also (e.g. JunD)?
# have a color bar beside plot marking whether each RS was statistically significant or not?

########
# make a ROC curve plot for EZH2/Suz12
