source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(curatedTCGAData)

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))

set.seed(1234)
# plotSubdir = "99_paper_figures/"
# scriptID = "99-paperFigures"
plotSubdir = "14_brca_DNAm_figures/"
createPlotSubdir(plotSubdir = plotSubdir)


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
# screen out region sets with less than 100 RS regions covered
keepInd = rsScores$regionSetCoverage >= 100

GRList = GRList[keepInd]
rsName = rsName[keepInd]
rsDescription = rsDescription[keepInd]
rsCollection = rsCollection[keepInd]
rsScores = rsScores[keepInd, ]

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
# Fig. 2a, PCA plot, DNA methylation

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

####################
# annoScoreDist


for (i in paste0("PC", 1:4)) {
    
    a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = i, 
                          pattern = c("esr1|eralpha", "gata3|foxa1|h3r17"), 
                          patternName = c("ER", "ER-related")) +
        theme(legend.position = c(0.15, 0.15)) +
        scale_color_manual(values = c("blue", "red", "orange")) + xlab(paste0("Region set rank (", i, ")"))
    a 
    ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", i, "_", .analysisID, ".svg"), 
           plot = a, device = "svg", width = plotWidth, height = plotHeight, units = plotUnits)
}

a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = "PC4", 
                      pattern = c("EZH2|SUZ12", "gata3|foxa1|h3r17"), 
                      patternName = c("Polycomb", "ER-related")) +
    theme(legend.position = c(0.15, 0.15)) +
    scale_color_manual(values = c("blue", "red", "orange")) + xlab("Region set rank (PC1)")
a

######################
# ER association by PC

corList = list()
thesePCs = paste0("PC", 1:4)
corVec = rep(-99, length(thesePCs))
names(corVec) = thesePCs
for (i in thesePCs) {
    # print(
    #     wilcox.test(pcaWithAnno[pcaWithAnno$ER_status == "Positive", i], 
    #                 pcaWithAnno[pcaWithAnno$ER_status == "Negative", i], conf.int = TRUE)
    # )
        corList[[i]]= cor.test(pcaWithAnno[, i], as.numeric(as.factor(pcaWithAnno$ER_status)) * 2 - 3, method = "spearman")
        corVec[i] = cor.test(pcaWithAnno[, i], 
                             as.numeric(as.factor(pcaWithAnno$ER_status)) * 2 - 3, 
                             method = "spearman")$estimate
}

corVec
medRank = rep(-99, length(thesePCs))
minRank = rep(-99, length(thesePCs))
for (i in seq_along(thesePCs)) {
    erInd = unique(c(grep(pattern = "esr1|eralpha", x = as.character(rsScores[order(rsScores[, thesePCs[i]], 
                                                        decreasing = TRUE), ]$rsName), ignore.case = TRUE),
      grep(pattern = "esr1|eralpha", x = as.character(rsScores[order(rsScores[, thesePCs[i]], 
                                                        decreasing = TRUE), ]$rsDescription), ignore.case = TRUE)))
    medRank[i]= median(erInd)
}

plot(abs(corVec), medRank)


####################
# interpretation of PC2
# EMT?

emtSig = as.character(read.csv(ffCode("COCOA_paper/metadata/EMT_core_sig_Groger2012.csv"), 
                  header = FALSE)[, 1])

rnaMAE = curatedTCGAData(diseaseCode = "BRCA", assays = c("RNASeq2GeneNorm"), 
                      dry.run = FALSE)
rna = assay(rnaMAE, "BRCA_RNASeq2GeneNorm-20160128")
colnames(rna) <- substr(colnames(rna), 1, 12)
sharedSamples = row.names(mPCA$x)[row.names(mPCA$x) %in% colnames(rna)]

rna <- rna[, sharedSamples]
mPCScores = mPCA$x[sharedSamples, ]

sGenes = intersect(emtSig, row.names(rna))
rna = rna[sGenes, ]

rnaPCA = prcomp(t(rna), center = TRUE, scale.=TRUE)
plot(rnaPCA$x[, c("PC1", "PC2")])
pcMean = colMeans(rnaPCA$x)
pcSD = apply(X = rnaPCA$x, 2, sd)
# zScorePC1 = abs((rnaPCA$x[, "PC1"] - pcMean[1]) / pcSD[1])
# zScorePC2 = abs((rnaPCA$x[, "PC2"] - pcMean[2]) / pcSD[2])
# sum(zScorePC1 > 3)
# sum(zScorePC2 > 3)
# outliers = row.names(rnaPCA$x[(zScorePC1 > 3) | (zScorePC2 > 3), ])

plot((rnaPCA$sdev^2 / sum(rnaPCA$sdev^2) )[1:10])

cor.test(rnaPCA$x[, "PC2"], mPCScores[, "PC2"])
all(row.names(mPCScores) == row.names(rnaPCA))


rnaD = dist(t(rna))
rnaClust = kmeans(x = t(rna), centers = 2)
table(rnaClust$cluster)
all(row.names(mPCScores) == names(rnaClust$cluster))
cor.test(x = mPCScores[, "PC2"], (rnaClust$cluster * 2 - 3), method = "spearman")
# not significant

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
                    signalCol=PCsToAnnotate, binNum=21, aggrMethod="default")

# for (i in seq_along(loadProfile)) {
#     ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "metaRegionPlots/", regionSetNames[i], "_", i), plot = profilePList[[i]], device = "pdf")
# }
ggsave(plot = mrProfileList$grob, 
       filename = ffPlot(paste0(plotSubdir, "metaRegionPlots", .analysisID, ".pdf")), device = "pdf")


# individual mr profiles

# normalize so plots from different PCs will be comparable
multiProfileP2 = normalizeMRProfile(signal=brcaCov, signalCol=signalCol, 
                       mrProfileList$metaRegionData, 
                       names(mrProfileList$metaRegionData))
multiProfileP2[["metaRegionData"]] = multiProfileP2
names(multiProfileP2[["metaRegionData"]])
topRSNames = c("wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak",
               "Human_MCF-7_H3R17me2_No-treatment_Brown.bed", 
               "wgEncodeAwgTfbsHaibT47dEraaV0416102Bpa1hUniPk.narrowPeak",
               "Human_MCF-7_FoxA1_No-treatment_Brown.bed",
               "GSM1501162_CEBPA.bed",
               "wgEncodeAwgTfbsSydhH1hescSuz12UcdUniPk.narrowPeak",
               "E104-H3K27me3.narrowPeak",
               "E032-H3K9me3.narrowPeak", 
               "wgEncodeAwgTfbsBroadH1hescEzh239875UniPk.narrowPeak")
abbrevNames = c("GATA3", "H3R17me2", "ER", "FOXA1", "CEBPA",
                "SUZ12", "H3K27me3", "H3K9me3",
                "EZH2")
# topRSNames = c("GSM835863_EP300.bed", 
#                "GSM607949_GATA1.bed")


for (i in seq_along(topRSNames)) {
    minVal = -1
    maxVal = 2
        
    pcP = multiProfileP2[["metaRegionData"]][topRSNames[i]]

    for (j in seq_along(signalCol[1:4])) {
        
        if (j == 4) {
            minVal = 0
            maxVal = 0.5
        } 
        
        
        myPlot = ggplot(data = filter(pcP[[1]], PC %in% signalCol[j]), mapping = aes(x =binID , y = loading_value)) + 
            # ggplot(data = pcP[[1]], mapping = aes(x =binID , y = loading_value)) + 
            geom_line() + ylim(c(minVal, maxVal)) + 
            # facet_wrap(facets = "PC") + 
            ggtitle(label = wrapper(topRSNames[i], width=30)) + xlab("Genome around Region Set, 14 kb") + 
            ylab("Normalized Correlation") + 
            theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
                  axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                  axis.title = element_blank(), title = element_blank(), axis.text.y=element_blank()
                  
            )
        myPlot
        ggsave(filename = ffPlot(paste0(plotSubdir, 
                                        "/mrProfilesWeightedMean_", abbrevNames[i],
                                        "_", signalCol[j], ".svg")), 
               plot = myPlot, device = "svg", height=20, width=30, units = "mm")
    }
    
    
}

################################################################################
# supplementary Fig
# DNA methylation in top region sets, ordered by PC score
# regions for which I plotted meta-region profiles

# get top region sets

for (i in seq_along(topRSNames)) {
    pdf(file = ffPlot(paste0(plotSubdir, "methylAlongPC1_", abbrevNames[i])))
    draw(signalAlongAxis(genomicSignal = brcaSharedC$methylProp, signalCoord = brcaCoord, 
                    regionSet = GRList[[topRSNames[i]]], 
                    sampleScores= mPCA$x[, c("PC1", "PC4")], orderByCol="PC1", 
                    topXVariables=100, variableScores = brcaCov[, "PC1"], 
                    cluster_columns=TRUE, show_row_names=FALSE)) 
    dev.off()
    pdf(file = ffPlot(paste0(plotSubdir, "methylAlongPC4_", abbrevNames[i])))
    draw(signalAlongAxis(genomicSignal = brcaSharedC$methylProp, signalCoord = brcaCoord, 
                         regionSet = GRList[[topRSNames[i]]], 
                         sampleScores= mPCA$x[, c("PC1", "PC4")], orderByCol="PC4", 
                         topXVariables=100, variableScores = brcaCov[, "PC4"], 
                         cluster_columns=TRUE, show_row_names=FALSE)) 
    dev.off()
}

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
