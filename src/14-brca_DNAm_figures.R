source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(curatedTCGAData)

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))

set.seed(1234)
# plotSubdir = "99_paper_figures/"
# scriptID = "99-paperFigures"
plotSubdir = "14_brca_DNAm_figures/"
createPlotSubdir(plotSubdir = plotSubdir)

# assign if not assigned so far
conAssign("dataID", "brcaDNAm657")
conAssign("nPerm", 300)
conAssign("variationMetric", "cov")

.analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)

plotWidth = 100
plotHeight = 100
plotUnits = "mm"

###########################################################
# reading in the region sets
# load LOLA database

loadGRList(genomeV = "hg38")


simpleCache(paste0("rsScores_", dataID, "_", variationMetric), 
            assignToVariable = "rsScores")
simpleCache(paste0("rsScores_", dataID, "_", variationMetric, "_median"), 
            assignToVariable="rsScoresMed")
simpleCache(paste0("rsScores_", dataID, "_loadings"), 
            assignToVariable = "rsScoresLoad")
# screen out region sets with less than 100 RS regions covered
keepInd = rsScores$regionSetCoverage >= 100

GRList = GRList[keepInd]
rsName = rsName[keepInd]
rsDescription = rsDescription[keepInd]
rsCollection = rsCollection[keepInd]
rsScores = rsScores[keepInd, ]
rsScoresMed = rsScoresMed[rsScoresMed$regionSetCoverage >= 100, ]
rsScoresLoad = rsScoresLoad[rsScoresLoad$regionSetCoverage >= 100, ]

#################################################################

allMPCAString = "allMPCA_657" #  "allMPCA_657"
# loading PCA and combining components that could separate ER=/-
# for rsEnrichment, PCs 1 and 4 could separate ER+/-
simpleCache(allMPCAString, assignToVariable = "mPCA")
brcaLoadings = mPCA$rotation


simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaSharedC", reload = TRUE)
brcaSharedC$methylProp = brcaSharedC$methylProp[, row.names(mPCA$x)]
brcaCoord = brcaSharedC[["coordinates"]]

# get top region sets for DNA methylation results
topScoreResultsInd = rsRankingIndex(rsScores, "PC1")$PC1[1:20]
topGRList = GRList[rsScores$rsName[topScoreResultsInd]]



# covariation for each CpG with PC
all(colnames(brcaSharedC$methylProp) == row.names(mPCA$x))
dim((mPCA$x))
dim(brcaSharedC$methylProp)
brcaCov = cov(t(brcaSharedC$methylProp), mPCA$x[, paste0("PC", 1:10)])

#############################################################################
# PCA plot, DNA methylation

colorByCols = "ER_status"
pcaWithAnno = cbind(mPCA$x, patientMetadata[row.names(mPCA$x) ,])

myPCList = list()
myPCList[[1]] = c("PC1", "PC2")
myPCList[[2]] = c("PC1", "PC3")
myPCList[[3]] = c("PC1", "PC4")
myPCList[[4]] = c("PC2", "PC3")
myPCList[[5]] = c("PC2", "PC4")
myPCList[[6]] = c("PC3", "PC4")

for (i in seq_along(myPCList)) {
    PCsToPlot = myPCList[[i]]
    multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                               plotCols = PCsToPlot, 
                                               colorByCols=colorByCols)
    ggplot2::ggsave(filename=ffPlot(paste0(plotSubdir,"/multiColorPCAPlots_allMPCA_", 
                                           PCsToPlot[1], "x", PCsToPlot[2], 
                                           ".pdf")), plot = multiColorPCAPlots, device = "pdf",
                    limitsize=FALSE)
}

# PC 1 and 2
aPCAPlot = ggplot(data = pcaWithAnno, mapping = aes(x = PC1, y= PC2)) + 
    geom_point(aes(col=ER_status), size = 1, alpha=0.3) + 
    theme(legend.position = c(0.15, 0.15), aspect.ratio = 1 
                          #axis.text = element_blank(), 
                          # axis.ticks = element_blank()
                          ) + scale_color_discrete(name="ER status") + 
    scale_x_continuous(breaks = scales::pretty_breaks(n=3)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n=3))
    
aPCAPlot
ggsave(filename = ffPlot(paste0(plotSubdir, "pc1_2_BRCA_DNAm.svg")), 
       height=(plotHeight/2)+10, width=(plotWidth/2)+10, units = "mm",
       plot = aPCAPlot, device = "svg")


# bigger one for presentation
biggerPlot = aPCAPlot = ggplot(data = pcaWithAnno, mapping = aes(x = PC1, y= PC2)) + 
    geom_point(aes(col=ER_status), size = 1, alpha=0.5) + 
    theme(legend.position = "right"
          #axis.text = element_blank(), 
          # axis.ticks = element_blank()
    ) + scale_color_discrete(name="ER status") + 
    scale_x_continuous(breaks = scales::pretty_breaks(n=3)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n=3)) + coord_equal()

biggerPlot

ggsave(filename = ffPlot(paste0(plotSubdir, "pc1_2_BRCA_DNAm_presentation.svg")), 
           height=80, width=100, units = "mm",
           plot = aPCAPlot, device = "svg")


pcaWithAnno = as.data.frame(pcaWithAnno)
for (i in c(paste0("PC", 1:4))) {
    print(
    wilcox.test(pcaWithAnno[pcaWithAnno$ER_status == "Positive", i], 
                pcaWithAnno[pcaWithAnno$ER_status == "Negative", i], conf.int = TRUE)
    )
    print(
    cor.test(pcaWithAnno[, i], as.numeric(as.factor(pcaWithAnno$ER_status)) * 2 - 3, method = "spearman")
    )
}

########
# get % variance explained for top PCs

varExpl = mPCA$sdev^2 / sum(mPCA$sdev^2) 
svg(filename = ffPlot(paste0(plotSubdir, "brca_DNAm_varianceExplained.svg")))
plot(varExpl)
dev.off()

plot(varExpl)
# 0.14010195 0.11371591 0.04713953 0.04632531 0.02099366 0.01939229

####################
# figure 2A
# annoScoreDist

for (i in paste0("PC", 1:4)) {
    
    a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = i, 
                          pattern = c("esr1|eralpha|eraa", "gata3|foxa1|h3r17", "ezh2|suz12"), 
                          patternName = c("ER", "ER-related", "polycomb")) +
        theme(legend.position = c(0.15, 0.15)) +
        scale_color_manual(values = c("blue", "red", "gray", "orange")) + 
        xlab(paste0("Region set rank (", i, ")")) + 
        theme(axis.title.y = element_blank(), 
              legend.text = element_blank(), 
              legend.title = element_blank(), 
              legend.position = "none") +
        scale_x_continuous(breaks = c(0, 1000, 2000), 
                         labels= c("0", "1000", "2000"), limits=c(-25, nrow(rsScores) + 25))
        
    a 
    ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", i, "_", .analysisID, ".svg"), 
           plot = a, device = "svg", width = plotWidth/2, height = plotHeight/2, units = plotUnits)
    
    ######################
    # using median scoring method
    a = plotAnnoScoreDist(rsScores = rsScoresMed, colsToPlot = i, 
                          pattern = c("esr1|eralpha|eraa", "gata3|foxa1|h3r17", "ezh2|suz12"), 
                          patternName = c("ER", "ER-related", "polycomb")) +
        theme(legend.position = c(0.15, 0.15)) +
        scale_color_manual(values = c("blue", "red", "gray", "orange")) + 
        xlab(paste0("Region set rank (", i, ")")) + 
        theme(axis.title.y = element_blank(), 
              legend.text = element_blank(), 
              legend.title = element_blank(), 
              legend.position = "none") +
        scale_x_continuous(breaks = c(0, 1000, 2000), 
                           labels= c("0", "1000", "2000"), limits=c(-25, nrow(rsScores) + 25))
    
    a 
    ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", i, "_", .analysisID, "_regionMedian", ".svg"), 
           plot = a, device = "svg", width = plotWidth/2, height = plotHeight/2, units = plotUnits)
    
}

# making a plot with legend
i="PC1"
a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = i, 
                      pattern = c("esr1|eralpha", "gata3|foxa1|h3r17", "ezh2|suz12"), 
                      patternName = c("ER", "ER-related", "Polycomb")) +
    theme(legend.position = c(0.15, 0.15)) +
    scale_color_manual(values = c("blue", "red", "gray", "orange")) + 
    xlab(paste0("Region set rank (", i, ")"))
    
    a 
ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", i, "_", .analysisID, "_withLegend.svg"), 
       plot = a, device = "svg", width = plotWidth, height = plotHeight, units = plotUnits)

# get correlation between median and mean scoring methods
# all(rsScores$rsName == rsScoresMed$rsName)
test = mapply(FUN = cor.test, x=rsScores[, paste0("PC", 1:4)], y = rsScoresMed[ , paste0("PC", 1:4)], method="spearman")
test["p.value", ] # p-value for PC1-4 is 0 (just 0 no decimalst)

test = mapply(FUN = cor.test, x=rsScoresLoad[, paste0("PC", 1:4)], y = rsScores[ , paste0("PC", 1:4)], method="pearson")
test["p.value", ] # p-value for PC1-4 is 0 (just 0 no decimalst)

###################
# rsScore distribution for COCOA poster
for (i in paste0("PC", 1:4)) {
    
    a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = i, 
                          pattern = c("esr1|eralpha|eraa", "gata3|foxa1|h3r17"), 
                          patternName = c("ER", "ER-related")) +
        theme(legend.position = c(0.15, 0.15)) +
        scale_color_manual(values = c("blue", "red", "gray")) + 
        xlab(paste0("Region set rank (", i, ")")) + 
        theme(axis.title.y = element_blank(), 
              legend.text = element_blank(), 
              legend.title = element_blank(), 
              legend.position = "none") +
        scale_x_continuous(breaks = c(0, 1000, 2000), 
                           labels= c("0", "1000", "2000"), limits=c(-25, nrow(rsScores) + 25))
    
    a 
    ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", i, "_Poster_", .analysisID, ".svg"), 
           plot = a, device = "svg", width = plotWidth/2, height = plotHeight/2, units = plotUnits)
}

i="PC1"
a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = i, 
                      pattern = c("esr1|eralpha", "gata3|foxa1|h3r17"), 
                      patternName = c("ER", "ER-related")) +
    theme(legend.position = c(0.15, 0.15)) +
    scale_color_manual(values = c("blue", "red", "gray")) + 
    xlab(paste0("Region set rank (", i, ")"))

a 
ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", i, "_Poster_", .analysisID, "_withLegend.svg"), 
       plot = a, device = "svg", width = plotWidth, height = plotHeight, units = plotUnits)



###################
# making a plot with just polycomb regions sets for supplement
for (i in paste0("PC", 1:10)) {
    
    a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = i, 
                          pattern = c("ezh2|suz12"), 
                          patternName = c("polycomb")) +
        theme(legend.position = c(0.15, 0.15)) +
        scale_color_manual(values = c("gray", "red")) + 
        xlab(paste0("Region set rank (", i, ")")) + 
        theme(axis.title.y = element_blank(), 
              legend.text = element_blank(), 
              legend.title = element_blank(), 
              legend.position = "none") +
        scale_x_continuous(breaks = c(0, 1000, 2000), 
                           labels= c("0", "1000", "2000"), limits=c(-25, nrow(rsScores) + 25))
    
    a 
    ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDistPolycomb_", i, "_", .analysisID, ".svg"), 
           plot = a, device = "svg", width = plotWidth/2, height = plotHeight/2, units = plotUnits)
}

a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = "PC4", 
                      pattern = c("EZH2|SUZ12", "gata3|foxa1|h3r17"), 
                      patternName = c("Polycomb", "ER-related")) +
    theme(legend.position = c(0.15, 0.15)) +
    scale_color_manual(values = c("blue", "red", "orange")) + xlab("Region set rank (PC1)")
a

######################
# figure 2B
# ER association by PC

corList = list()
thesePCs = paste0("PC", 1:4)
corVec = rep(-99, length(thesePCs))
names(corVec) = thesePCs
tmp = rep(-99, length(thesePCs))
resultsDF = data.frame(spearmanCor=tmp,
                       spearmanPVal=tmp,
                       medianDiff=tmp, 
                       wilcoxPVal=tmp)
row.names(resultsDF) = thesePCs

for (i in thesePCs) {
    
        tmp = wilcox.test(pcaWithAnno[pcaWithAnno$ER_status == "Positive", i],
                    pcaWithAnno[pcaWithAnno$ER_status == "Negative", i], conf.int = TRUE)
        resultsDF[i, "medianDiff"] = tmp$estimate
        resultsDF[i, "wilcoxPVal"] = tmp$p.value
        
        
        corList[[i]]= cor.test(pcaWithAnno[, i], as.numeric(as.factor(pcaWithAnno$ER_status)) * 2 - 3, method = "spearman")
        corVec[i] = cor.test(pcaWithAnno[, i], 
                             as.numeric(as.factor(pcaWithAnno$ER_status)) * 2 - 3, 
                             method = "spearman")$estimate
        tmp2 =  cor.test(pcaWithAnno[, i], 
                         as.numeric(as.factor(pcaWithAnno$ER_status)) * 2 - 3, 
                         method = "spearman")
        resultsDF[i, "spearmanCor"] = tmp2$estimate
        resultsDF[i, "spearmanPVal"] = tmp2$p.value
}

pcSD = apply(X = pcaWithAnno[, thesePCs], MARGIN = 2, FUN = sd)
resultsDF$normMedianDiff = resultsDF$medianDiff / pcSD
resultsDF = cbind(PC=thesePCs, resultsDF)
resultsDF
corVec
write.csv(x = resultsDF, file = ffSheets(paste0("pcERStatusRelationship","_", 
                                                         dataID, ".csv")),row.names = FALSE)
tmedRank = rep(-99, length(thesePCs))
minRank = rep(-99, length(thesePCs))
for (i in seq_along(thesePCs)) {
    erInd = unique(c(grep(pattern = "esr1|eralpha", x = as.character(rsScores[order(rsScores[, thesePCs[i]], 
                                                        decreasing = TRUE), ]$rsName), ignore.case = TRUE),
      grep(pattern = "esr1|eralpha", x = as.character(rsScores[order(rsScores[, thesePCs[i]], 
                                                        decreasing = TRUE), ]$rsDescription), ignore.case = TRUE)))
    medRank[i]= median(erInd)
}

plot(abs(corVec), medRank)

###
# order samples by PC score, color by ER status
erStatusDF = cbind(apply(X = mPCA$x[, paste0("PC", 1:4)], 
                         MARGIN = 2, FUN = function(x) order(order(x, decreasing = FALSE))), 
                   as.data.frame(patientMetadata[row.names(mPCA$x) , ])[, c("subject_ID", "ER_status")])
erStatusDF = pivot_longer(erStatusDF, cols = c("PC1", "PC2", "PC3", "PC4"), names_to = "PC", values_to = "rank")    

erStatusDF$barHeight = rep(1, nrow(erStatusDF))
erStatusDF = arrange(erStatusDF, PC, rank) 

erStatusPlot = ggplot(data = erStatusDF, mapping = aes(x=rank, y=barHeight, group=PC)) + 
    geom_col(aes(col=ER_status)) + scale_color_discrete(breaks=c("Positive", "Negative")) +
    xlab("Samples ordered by PC score") + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                             axis.title.y = element_blank(), axis.line = element_blank())
erStatusPlot

ggplot2::ggsave(filename=ffPlot(paste0(plotSubdir,"/orderedERStatus.svg")), 
                plot = erStatusPlot, device = "svg", height = plotHeight / 2, 
                width = plotWidth, units = plotUnits)


####################
# # interpretation of PC2
# # EMT?
# 
# emtSig = as.character(read.csv(ffCode("COCOA_paper/metadata/EMT_core_sig_Groger2012.csv"), 
#                   header = FALSE)[, 1])
# 
# rnaMAE = curatedTCGAData(diseaseCode = "BRCA", assays = c("RNASeq2GeneNorm"), 
#                       dry.run = FALSE)
# rna = assay(rnaMAE, "BRCA_RNASeq2GeneNorm-20160128")
# colnames(rna) <- substr(colnames(rna), 1, 12)
# sharedSamples = row.names(mPCA$x)[row.names(mPCA$x) %in% colnames(rna)]
# 
# rna <- rna[, sharedSamples]
# mPCScores = mPCA$x[sharedSamples, ]
# 
# sGenes = intersect(emtSig, row.names(rna))
# rna = rna[sGenes, ]
# 
# rnaPCA = prcomp(t(rna), center = TRUE, scale.=TRUE)
# plot(rnaPCA$x[, c("PC1", "PC2")])
# pcMean = colMeans(rnaPCA$x)
# pcSD = apply(X = rnaPCA$x, 2, sd)
# # zScorePC1 = abs((rnaPCA$x[, "PC1"] - pcMean[1]) / pcSD[1])
# # zScorePC2 = abs((rnaPCA$x[, "PC2"] - pcMean[2]) / pcSD[2])
# # sum(zScorePC1 > 3)
# # sum(zScorePC2 > 3)
# # outliers = row.names(rnaPCA$x[(zScorePC1 > 3) | (zScorePC2 > 3), ])
# 
# plot((rnaPCA$sdev^2 / sum(rnaPCA$sdev^2) )[1:10])
# 
# cor.test(rnaPCA$x[, "PC2"], mPCScores[, "PC2"])
# all(row.names(mPCScores) == row.names(rnaPCA))
# 
# 
# rnaD = dist(t(rna))
# rnaClust = kmeans(x = t(rna), centers = 2)
# table(rnaClust$cluster)
# all(row.names(mPCScores) == names(rnaClust$cluster))
# cor.test(x = mPCScores[, "PC2"], (rnaClust$cluster * 2 - 3), method = "spearman")
# # not significant

###################
# Fig. 2c, meta region loading profile plots, DNA methylation BRCA
# Supplementary fig. 

PCsToAnnotate = paste0("PC", 1:10)
signalCol = PCsToAnnotate
# # manually get region set name
# regionSetNames = c("Gata3", "H3R17me2", "ESR1", "Gata3", 
#                    "FoxA1", "ESR1", "ESR1", "ESR1", 
#                    "FoxA1", "ESR1", "AR", "AR", 
#                    "Znf217", "TCF7L2", "ESR1", "ESR1", 
#                    "JunD", "ESR1", "FoxA1", "H3R17me2")
# get top region sets for DNA methylation results
topScoreResultsInd = rsRankingIndex(rsScores, "PC1")$PC1[1:50]
topScoreResultsInd = c(topScoreResultsInd, rsRankingIndex(rsScores, "PC4")$PC4[1:20])
topGRList = GRList[rsScores$rsName[topScoreResultsInd]]

regionSetNames = names(topGRList)
wideGRList <- lapply(topGRList, resize, width=14000, fix="center")

mrProfileList  = makeMetaRegionPlots(signal=brcaCov, signalCoord=brcaCoord,
                    GRList=wideGRList, rsNames=regionSetNames, 
                    signalCol=PCsToAnnotate, binNum=21, aggrMethod="default", 
                    normMethod = "none")

# for (i in seq_along(loadProfile)) {
#     ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "metaRegionPlots/", regionSetNames[i], "_", i), plot = profilePList[[i]], device = "pdf")
# }
ggsave(plot = mrProfileList$grob, 
       filename = ffPlot(paste0(plotSubdir, "metaRegionPlots", .analysisID, ".pdf")), device = "pdf")


# individual mr profiles

# normalize so plots from different PCs will be comparable
# multiProfileP2 = normalizeMRProfile(signal=brcaCov, signalCol=signalCol, 
#                        mrProfileList$metaRegionData, 
#                        names(mrProfileList$metaRegionData), 
#                        normMethod = "mean")
multiProfileP2 = normalizeMRProfile(signal=brcaCov, signalCol=signalCol, 
                                    mrProfileList$metaRegionData, 
                                    names(mrProfileList$metaRegionData),
                                    normMethod = "zscore")
# multiProfileP2 = normalizeMRProfile(signal=brcaCov, signalCol=signalCol, 
#                                     mrProfileList$metaRegionData, 
#                                     names(mrProfileList$metaRegionData),
#                                     normMethod = "normPVal")

multiProfileP2[["metaRegionData"]] = multiProfileP2
names(multiProfileP2[["metaRegionData"]])
topRSNames = c("wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak",
               "Human_MCF-7_H3R17me2_No-treatment_Brown.bed", 
               "wgEncodeAwgTfbsHaibT47dEraaV0416102Bpa1hUniPk.narrowPeak",
               "Human_MCF-7_FoxA1_No-treatment_Brown.bed",
               "GSM1501162_CEBPA.bed",
               "GSM614003_TAL1.bed",
               "wgEncodeAwgTfbsSydhH1hescSuz12UcdUniPk.narrowPeak",
               "E104-H3K27me3.narrowPeak",
               "E032-H3K9me3.narrowPeak", 
               "wgEncodeAwgTfbsBroadH1hescEzh239875UniPk.narrowPeak")
abbrevNames = c("GATA3", "H3R17me2", "ER", "FOXA1", 
                "SUZ12", "H3K27me3", "H3K9me3",
                "EZH2")
# topRSNames = c("GSM835863_EP300.bed", 
#                "GSM607949_GATA1.bed")

for (i in seq_along(topRSNames)) {
    minVal = -0.5
    maxVal = 2.5
        
    pcP = multiProfileP2[["metaRegionData"]][topRSNames[i]]

    for (j in seq_along(signalCol[1:4])) {
        
        # if (j == 4) {
        #     minVal = 0
        #     maxVal = 0.5
        # } 
        
        myPlot = ggplot(data = filter(pcP[[1]], PC %in% signalCol[j]), mapping = aes(x =binID , y = loading_value)) + 
            # ggplot(data = pcP[[1]], mapping = aes(x =binID , y = loading_value)) + 
            geom_line() + ylim(c(minVal, maxVal)) + geom_hline(yintercept = 0, col="red", alpha = 0.25) +
            # facet_wrap(facets = "PC") + 
            ggtitle(label = wrapper(topRSNames[i], width=30)) + xlab("Genome around Region Set, 14 kb") + 
            ylab("Normalized Correlation") + 
            theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
                  axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                  axis.title = element_blank(), title = element_blank(), axis.text.y=element_blank()) 
        myPlot
        ggsave(filename = ffPlot(paste0(plotSubdir, 
                                        "/mrProfiles_", abbrevNames[i],
                                        "_", signalCol[j], ".svg")), 
               plot = myPlot, device = "svg", height=20, width=30, units = "mm")
    }
}
################################################################################
# metaregion plots for the median score

mrProfileList  = makeMetaRegionPlots(signal=brcaCov, signalCoord=brcaCoord,
                                     GRList=wideGRList, rsNames=regionSetNames, 
                                     signalCol=PCsToAnnotate, binNum=21, aggrMethod="regionMedian", 
                                     normMethod = "none")

# for (i in seq_along(loadProfile)) {
#     ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "metaRegionPlots/", regionSetNames[i], "_", i), plot = profilePList[[i]], device = "pdf")
# }
ggsave(plot = mrProfileList$grob, 
       filename = ffPlot(paste0(plotSubdir, "metaRegionPlots", .analysisID, "_regionMedian.pdf")), device = "pdf")


# individual mr profiles

# normalize so plots from different PCs will be comparable
# multiProfileP2 = normalizeMRProfile(signal=brcaCov, signalCol=signalCol, 
#                        mrProfileList$metaRegionData, 
#                        names(mrProfileList$metaRegionData), 
#                        normMethod = "mean")
multiProfileP2 = normalizeMRProfile(signal=brcaCov, signalCol=signalCol, 
                                    mrProfileList$metaRegionData, 
                                    names(mrProfileList$metaRegionData),
                                    normMethod = "medianZScore")
# multiProfileP2 = normalizeMRProfile(signal=brcaCov, signalCol=signalCol, 
#                                     mrProfileList$metaRegionData, 
#                                     names(mrProfileList$metaRegionData),
#                                     normMethod = "normPVal")

multiProfileP2[["metaRegionData"]] = multiProfileP2
names(multiProfileP2[["metaRegionData"]])
topRSNames = c("wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak",
               "Human_MCF-7_H3R17me2_No-treatment_Brown.bed", 
               "wgEncodeAwgTfbsHaibT47dEraaV0416102Bpa1hUniPk.narrowPeak",
               "Human_MCF-7_FoxA1_No-treatment_Brown.bed",
               "GSM1501162_CEBPA.bed",
               "GSM614003_TAL1.bed",
               "wgEncodeAwgTfbsSydhH1hescSuz12UcdUniPk.narrowPeak",
               "E104-H3K27me3.narrowPeak",
               "E032-H3K9me3.narrowPeak", 
               "wgEncodeAwgTfbsBroadH1hescEzh239875UniPk.narrowPeak")
abbrevNames = c("GATA3", "H3R17me2", "ER", "FOXA1", 
                "SUZ12", "H3K27me3", "H3K9me3",
                "EZH2")
# topRSNames = c("GSM835863_EP300.bed", 
#                "GSM607949_GATA1.bed")

for (i in seq_along(topRSNames)) {
    minVal = -0.5
    maxVal = 4.5
    
    pcP = multiProfileP2[["metaRegionData"]][topRSNames[i]]
    
    for (j in seq_along(signalCol[1:4])) {
        
        # if (j == 4) {
        #     minVal = 0
        #     maxVal = 0.5
        # } 
        
        myPlot = ggplot(data = filter(pcP[[1]], PC %in% signalCol[j]), mapping = aes(x =binID , y = loading_value)) + 
            # ggplot(data = pcP[[1]], mapping = aes(x =binID , y = loading_value)) + 
            geom_line() + ylim(c(minVal, maxVal)) + geom_hline(yintercept = 0, col="red", alpha = 0.25) +
            # facet_wrap(facets = "PC") + 
            ggtitle(label = wrapper(topRSNames[i], width=30)) + xlab("Genome around Region Set, 14 kb") + 
            ylab("Normalized Covariation") + 
            theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
                  axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                  axis.title = element_blank(), title = element_blank(), axis.text.y=element_blank()) 
        myPlot
        ggsave(filename = ffPlot(paste0(plotSubdir, 
                                        "/mrProfiles_", abbrevNames[i],
                                        "_", signalCol[j], "_regionMedian.svg")), 
               plot = myPlot, device = "svg", height=20, width=30, units = "mm")
    }
}

###############################################################################
# get p-values 
# this function currently requires that an odd number of bins was used (middle
# bin is selected)
# @param binDF one row per bin. Columns: average, coverage, standard deviation
# @param ? test type?
getMRPVal <- function(binDF, ) {
    
    if (nrow(binDF)) {
        # compare means
        
        # need n (coverage) and standard deviation
        
        # wilcoxon rank sum test
    }

}

####################
# supplementary figure
# histone modification region sets: e.g. H3K9me3, H3K27me3

# top 2 H3K9me3 and H3K27me3 from PC2
topRSNames = c("E122-H3K9me3.narrowPeak",
               "E021-H3K9me3.narrowPeak",
               "E117-H3K27me3.narrowPeak",
               "E104-H3K27me3.narrowPeak",
               "GSM1501162_CEBPA.bed",
               "GSM614003_TAL1.bed")
               
abbrevNames = c("HUVEC H3K9me3", "iPS_DF_6.9_Cell_Line H3K9me3", "HeLa_H3K27me3",
                "Right Atrium H3K27me3", "CEBPA", "TAL1")

topGRList = GRList[topRSNames]

regionSetNames = names(topGRList)
wideGRList <- lapply(topGRList, resize, width=14000, fix="center")

mrProfileList  = makeMetaRegionPlots(signal=brcaCov, signalCoord=brcaCoord,
                                     GRList=wideGRList, rsNames=regionSetNames, 
                                     signalCol=PCsToAnnotate, binNum=21, aggrMethod="default", normMethod = "none")

multiProfileP2 = normalizeMRProfile(signal=brcaCov, signalCol=signalCol, 
                                    mrProfileList$metaRegionData, 
                                    names(mrProfileList$metaRegionData),
                                    normMethod = "zscore")
multiProfileP2[["metaRegionData"]] = multiProfileP2

signalCol = paste0("PC", 1:4)
for (i in seq_along(topRSNames)) {
    minVal = -0.5
    maxVal = 2.5
    
    pcP = multiProfileP2[["metaRegionData"]][topRSNames[i]]
    
    myPlot = ggplot(data = filter(pcP[[1]], PC %in% signalCol), mapping = aes(x =binID , y = loading_value)) + 
        theme(panel.border = element_rect(fill = NA)) +
        geom_line() + ylim(c(minVal, maxVal)) + geom_hline(yintercept = 0, col="red", alpha = 0.25) +
        facet_wrap(facets = "PC", nrow = length(signalCol), ncol = 1) + 
        ggtitle(label = wrapper(topRSNames[i], width=30)) + xlab("Genome around Region Set, 14 kb") + 
        ylab("Normalized Covariance") + 
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
              axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.title = element_blank(), title = element_blank(), axis.text.y=element_blank(),
              strip.background = element_blank(), strip.text = element_blank()) 
    myPlot
    ggsave(filename = ffPlot(paste0(plotSubdir, 
                                    "/mrProfilesWeightedMean_", abbrevNames[i],
                                    ".svg")), 
           plot = myPlot, device = "svg", height=100, width=25, units = "mm")
}

# one plot with legend
myPlot = myPlot + theme(axis.text.y = element_text())
ggsave(filename = ffPlot(paste0(plotSubdir, 
                                "/mrProfilesWeightedMean_", abbrevNames[i],
                                "_withLegend.svg")), 
       plot = myPlot, device = "svg", height=100, width=25, units = "mm")



################################################################################
# supplementary Fig
# DNA methylation in top region sets, ordered by PC score
# regions for which I plotted meta-region profiles

# get top region sets
colnames(brcaSharedC$methylProp) = gsub(pattern = "-",
                                        replacement = "_", 
                                        x = colnames(brcaSharedC$methylProp))
row.names(mPCA$x) = gsub(pattern = "-",
                        replacement = "_", 
                        x = row.names(mPCA$x))
myPCs = paste0("PC", 1:4)
for (i in seq_along(topRSNames)) {
    thisRSM = COCOA:::averagePerRegion(signal=brcaSharedC$methylProp, 
                                       signalCoord = brcaCoord, 
                                       regionSet = GRList[[topRSNames[i]]], 
                                       signalCol = colnames(brcaSharedC$methylProp))
    thisRSCovScores = COCOA:::averagePerRegion(signal=brcaCov[, myPCs], 
                                               signalCoord = brcaCoord, 
                                               regionSet = GRList[[topRSNames[i]]], 
                                               signalCol = myPCs)
    thisRSM = as.data.frame(thisRSM)
    thisRSCovScores = as.data.frame(thisRSCovScores)
    for (j in seq_along(myPCs)) {
        
        # pdf(file = ffPlot(paste0(plotSubdir, 
        #                          "methylAlong", myPCs[j], "_", abbrevNames[i], ".pdf")))
        png(filename = ffPlot(paste0(plotSubdir,
                                     "methylAlong", myPCs[j], "_", abbrevNames[i], ".png")),
            width = 100, height = 75, units = "mm", pointsize = 12,
            bg = "white",  res = 300,
            type = c("cairo", "cairo-png", "Xlib", "quartz"))
        # # svg(filename = ffPlot(paste0(plotSubdir, "methylAlong", myPCs[j], "_", abbrevNames[i], ".svg")))
        draw(signalAlongAxis(genomicSignal = thisRSM[, !(colnames(thisRSM) %in% c("chr", "start", "end"))], 
                             signalCoord = thisRSM[, c("chr", "start", "end")], 
                             regionSet = GRList[[topRSNames[i]]], 
                             sampleScores= mPCA$x[, myPCs], orderByCol=myPCs[j], 
                             topXVariables=100, variableScores = thisRSCovScores[, myPCs[j]], 
                             cluster_columns=TRUE, show_row_names=FALSE, show_column_names=FALSE, 
                             column_title = "Regions (top 100)", name = "DNA methylation level", 
                             row_title = "Samples ordered by PC score", 
                             row_title_gp = gpar(fontsize = 14), # 54
                             column_title_gp = gpar(fontsize = 14),
                             show_heatmap_legend = FALSE
                             )
             ) 
        dev.off()
        
        if ((j == 1) && (i == 1)) {
            svg(filename = ffPlot(paste0(plotSubdir, "methylAlong", myPCs[j], "_", abbrevNames[i], "_hasLegend.svg")), 
                width = 10, height=10) # inches
            draw(signalAlongAxis(genomicSignal = thisRSM[, !(colnames(thisRSM) %in% c("chr", "start", "end"))], 
                                 signalCoord = thisRSM[, c("chr", "start", "end")], 
                                 regionSet = GRList[[topRSNames[i]]], 
                                 sampleScores= mPCA$x[, myPCs], orderByCol=myPCs[j], 
                                 topXVariables=10, variableScores = thisRSCovScores[, myPCs[j]], 
                                 cluster_columns=TRUE, show_row_names=FALSE, show_column_names=FALSE, 
                                 column_title = "Regions (top 100)", name = "DNA methylation level", 
                                 row_title = "Samples ordered by PC score")) 
            dev.off()
            # row_title_gp = gpar(fontsize = 14), # 54
            # column_title_gp = gpar(fontsize = 14))
        }

    }
}

# make one plot as svg, with the legend 

########################################################################################
# for supplementary fig comparing mean to median scoring
# also for supplementary fig comparing scoring with loadings to using covariance values

signalCol = paste0("PC", 1:10)

# convert to long format
lRSScores = tidyr::pivot_longer(rsScores, cols = all_of(signalCol), names_to = "PC", values_to = "meanScore")
# lRSScores$scoringMetric = "mean"

lRSScoresMed = tidyr::pivot_longer(rsScoresMed, cols = all_of(signalCol), names_to = "PC", values_to = "medianScore")
# lRSScoresMed$scoringMetric = "median"
lRSScoresLoad =  tidyr::pivot_longer(rsScoresLoad, cols = all_of(signalCol), names_to = "PC", values_to = "loadings")

all(lRSScores$rsName == lRSScoresMed$rsName)


lRSScoresBoth = inner_join(lRSScores, lRSScoresMed, by=c("rsName", "rsDescription", "PC"))
lRSScoresBoth = inner_join(lRSScoresBoth, lRSScoresLoad, by=c("rsName", "rsDescription", "PC"))

# annotation for color in plots
lRSScoresBoth$rsGroup = NA


myPattern = c("esr1|eralpha|eraa", "gata3|foxa1|h3r17", "ezh2|suz12")
patternName = c("ER", "ER-related", "polycomb")
for (i in seq_along(myPattern)) {
    lRSScoresBoth$rsGroup[grep(pattern = myPattern[i], x = lRSScoresBoth$rsName, 
                               ignore.case = TRUE)] = patternName[i]
    lRSScoresBoth$rsGroup[grep(pattern = myPattern[i], x = lRSScoresBoth$rsDescription, 
                               ignore.case = TRUE)] = patternName[i]
                                                                      
}
lRSScoresBoth$rsGroup = as.factor(lRSScoresBoth$rsGroup)

# color by region set category

for (i in 1:4) {
    a = ggplot(data = filter(lRSScoresBoth, PC == paste0("PC", i)), mapping = aes(x=meanScore, y=medianScore)) + 
        geom_point(aes(col=rsGroup), alpha=0.5, shape=3) +
        scale_color_manual(values = c("blue", "red", "orange"), na.value=alpha("gray", 0.01)) + 
        xlab("Mean scoring method") + ylab("Median scoring method") +
        theme(legend.text = element_blank(),
              legend.title = element_blank(),
              legend.position = "none")
    
    ggplot2::ggsave(filename=ffPlot(paste0(plotSubdir,"/median_mean_correlation_", signalCol[i], 
                                           ".svg")), plot = a, device = "svg",
                    width =  plotWidth/2, height = plotHeight/2, units = plotUnits)
    
    a = ggplot(data = filter(lRSScoresBoth, PC == paste0("PC", i)), mapping = aes(x=meanScore, y=loadings)) + 
        geom_point(aes(col=rsGroup), alpha=0.75, shape=3) +
        scale_color_manual(values = c("blue", "red", "orange"), na.value=alpha("gray", alpha = 0.0001)) + 
        xlab("Using covariance") + ylab("Using loadings") +
        theme(legend.text = element_blank(),
              legend.title = element_blank(),
              legend.position = "none")
    
    ggplot2::ggsave(filename=ffPlot(paste0(plotSubdir,"/cor_loadings_correlation_", signalCol[i], 
                                           ".svg")), plot = a, device = "svg",
                    width =  plotWidth/1.5, height = plotHeight/1.5, units = plotUnits)
    

    
}

# legend

a = ggplot(data = filter(lRSScoresBoth, PC == paste0("PC", i)), mapping = aes(x=meanScore, y=medianScore)) + 
    geom_point(aes(col=rsGroup), alpha=0.5, shape=3) +
    theme(legend.position = c(0.8, 0.1)) +
    scale_color_manual(values = c("blue", "red", "orange"), na.value="gray") + 
    xlab("Mean scoring method") + ylab("Median scoring method")

ggplot2::ggsave(filename=ffPlot(paste0(plotSubdir,"/median_mean_correlation_", signalCol[i], 
                                       "_withLegend.svg")), plot = a, device = "svg")
t########################################################################################
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

#########################################################################################
# exploratory visualization

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


###############################

a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = "PC1", pattern = "esr1|eralpha")
a + scale_color_brewer(palette = "Dark2")
a + scale_color_brewer(palette = "RdYlBu")
a + scale_color_brewer(palette = "Paired")
a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = "PC1", pattern = c("gata3|foxa1|H3R17me2", "esr1|eralpha"))
a + scale_color_brewer(palette = "Dark2")
a + scale_color_brewer(palette = "RdYlBu")
a + scale_color_brewer(palette = "Paired")
