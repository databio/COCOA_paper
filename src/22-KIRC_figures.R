
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(cowplot)
library(ROCR)

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


plotAnnoScoreDist <- function(rsScores, colsToPlot, pattern, patternName=pattern) {
    
    realRSScores = rsScores
    
    realRSScores$rank = order(order(as.numeric(realRSScores[, colsToPlot]), decreasing = TRUE))
    
    rsGroup = rep(0, nrow(realRSScores))
    for (i in seq_along(pattern)) {
        thisPatternInd = grepl(pattern = pattern[i], 
                               x = realRSScores$rsName, 
                               ignore.case = TRUE)
        # earlier groups will be overwritten if they are not mutually exclusive
        # (e.g. if the same region set matches two patterns, the later group
        # will be assigned)
        rsGroup[thisPatternInd] = i
    }
    rsGroup = as.factor(rsGroup)
    # necessary so region set names will be in legend
    levels(rsGroup) = c("Other", patternName)
    realRSScores$Group = rsGroup 
    
    
    rsCorP = ggplot(data = realRSScores, 
                    mapping = aes(x=rank, y=get(colsToPlot), 
                                  col=Group), alpha=0.1) + # alpha=Group
        geom_point() +
        ylab("Region set score") + xlab("Region set rank")
    # add each group (Other and pattern) sequentially so they will be plotted on top
    # of other points ('Other' plotted first)

    for (i in 2:(length(pattern) + 1)) {
        rsCorP = rsCorP + geom_point(data = realRSScores[as.numeric(realRSScores$Group) == i, ])
    }
    
    # rsCorP + coord_flip()
    polycombStatusP = ggplot(data = realRSScores, mapping = aes(x=rank, y=as.numeric(Group)-1)) + 
        geom_col(aes(col=Group)) 
    polycombStatusP 
    
    gg_dist_g1 = polycombStatusP
    gg_dist_g2 = rsCorP
    gg_scatter = rsCorP
    
    # gg_dist_g2 = gg_dist_g2 + coord_flip()
    
    # Remove some duplicate axes
    gg_dist_g1 = gg_dist_g1 + theme(axis.title.x=element_blank(),
                                    axis.title.y = element_blank(),
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
    p = first_col
    p
    
    return(p)
}

a = plotAnnoScoreDist(rsScores = realRSScores, colsToPlot = "cancerStage", pattern = c("EZH2|SUZ12", "H3K27me3"))
a + ggtitle("Association between DNA methylation and cancer stage ")
plotAnnoScoreDist(rsScores = realRSScores, colsToPlot = "cancerStage", pattern = c("EZH2|SUZ12"))

# either color dots for EZH2/Suz12 or have color bar beside plot
# make region sets that I will be pulling out to plot separately a different color also (e.g. JunD)?
# have a color bar beside plot marking whether each RS was statistically significant or not?

########
# make a ROC curve plot for EZH2/Suz12
pred = realRSScores$cancerStage / max(realRSScores$cancerStage, na.rm = TRUE)
rocPreds = ROCR::prediction(pred, grepl(pattern = "EZH2|SUZ12", 
                                        x = realRSScores$rsName, 
                                        ignore.case = TRUE))
testAUC = ROCR::performance(rocPreds, measure="auc")@y.values[[1]]
testAUC
perf = ROCR::performance(rocPreds, measure = "tpr", x.measure = "fpr")
plot(perf)
title(main= "ROC curve for EZH2 and SUZ12")
