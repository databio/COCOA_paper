
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(MIRA)
library(MultiAssayExperiment)
library(ensembldb)
library(EnsDb.Hsapiens.v75)

if (!exists("dataID")) {
    dataID = "CLL196MOFA"    
}
if (!exists("variationMetric")) {
    variationMetric = "cov"   
}
plotSubdir = "33-MOFA_figs/"
if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

inputID = dataID

##############################################################################
# Fig. 4, MOFA CLL
# http://msb.embopress.org/content/14/6/e8124
# factors 1 (strong), 7 and 9 were associated with DNA methylation
# factor 1: cell type/differentiation, factor 7: chemo-immunotherapy treatment prior to sample collection
# factor 7: del17p, TP53 mutations, methylation of oncogenes


loadMOFAData(methylMat = TRUE, signalCoord=TRUE, latentFactorMat = TRUE,
             cllMultiOmics=TRUE, featureLFCor = TRUE)
simpleCache("mofaPermZScoresCor", assignToVariable = "rsZScores")
simpleCache(paste0("rsScore_", dataID, "_", variationMetric), assignToVariable = "rsScores")


# get same order
MAE = MultiAssayExperiment(list(mutations = cllMultiOmics$Mutation, 
                                latentFactors=t(latentFactorMat)))
mutDF = as.data.frame(assay(MAE[, complete.cases(MAE), ], "mutations"))
latentFactors = t(assay(MAE[, complete.cases(MAE), ], "latentFactors"))
IGHV = mutDF["IGHV", ]
IGHV[IGHV == 1] = "mutated"
IGHV[IGHV == 0] = "unmutated"
IGHV = as.factor(IGHV)
latentFactors = data.frame(latentFactors, IGHV, KLHL6=mutDF["KLHL6", ], BRAF=factor(mutDF["BRAF", ])) 
# checking whether LF9 is associated with certain immunce cell markers
# CD56=cllMultiOmics$mRNA["ENSG00000149294", ],
# CD3=cllMultiOmics$mRNA["ENSG00000167286", ])
# CD19=cllMultiOmics$mRNA["ENSG00000177455", ]
# latentFactors$IGHV = factor(latentFactors$IGHV, levels = c("mutated", "unmutated"))
# cor.test(x = latentFactors$CD3, y = latentFactors$LF9)

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


#################### plot raw data in top regions for LF1
# simpleCache("inferredMethylWeightsMOFA", assignToVariable = "featureLFCor")
# make sure they are in the same order/have same CpGs
featureLFCor = featureLFCor[row.names(methylMat), ]
loadGRList(genomeV = "hg19")
dataID = "cll196"
lfCols= paste0("LF", 1:10)
# lfCols= paste0("LF", c(1:3, 5:7, 9))
cpgCov = rsScores$signalCoverage
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
        print(signalAlongAxis(genomicSignal=methylMat,
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

signalAlongAxis(genomicSignal=methylMat,
              signalCoord=signalCoord,
              sampleScores=latentFactorMat,
              regionSet=topLF1RS[["Gm12878_WE.bed"]], orderByCol=lfCols[i],
              topXVariables=50,
              variableScores = abs(as.numeric(featureLFCor[, lfCols[i]])),
              cluster_columns = TRUE, column_title = "Individual cytosines")

###############################################################################
# Panel B
# sidebar with IGHV status for methylation heatmap
latentFactors = latentFactors[order(latentFactors$LF1, decreasing = TRUE), ]
latentFactors$LF1Index = 1:nrow(latentFactors)
latentFactors$IGHV = factor(latentFactors$IGHV, levels = c("unmutated", "mutated"))

IGHVBar = ggplot(data = latentFactors, mapping = aes(x = LF1Index, y = 1)) + 
    geom_col(aes(color=IGHV)) + coord_flip() + scale_fill_discrete(na.value = "gray") + 
    theme(axis.title = element_blank(), axis.text = element_blank(), 
          axis.line = element_blank(), axis.ticks = element_blank())

IGHVBar
ggsave(filename = ffPlot(paste0(plotSubdir, "IGHVBarLF1.svg")), 
       plot = IGHVBar, 
       device = "svg", width = 50, height = 35, units = "mm")

##############################################################################
# figure connecting LF2 and chr12 trisomy to NANOG
# NANOG is ENSG00000111704
library(BloodCancerMultiOmics2017)
loadMOFAData(methylMat = TRUE, signalCoord=TRUE, latentFactorMat = TRUE,
             cllMultiOmics=TRUE, featureLFCor = TRUE)
data("dds", package = "BloodCancerMultiOmics2017")
cllRNA = assay(dds)
colnames(cllRNA) <- colData(dds)$PatID
head(cllRNA)

# thisRegionSet = "GSM1124068_SOX2.bed"
# thisGeneInd = grep(pattern = "ENSG00000181449", x = row.names(cllRNA))
thisRegionSet = "GSM1124071_NANOG.bed"
thisGeneInd = grep(pattern = "ENSG00000111704", x = row.names(cllRNA))

hist(cllRNA[thisGeneInd, ])
rnaIDs = colnames(cllRNA)

# match ID's 
sharedIDs = rnaIDs[rnaIDs %in% row.names(latentFactorMat)]

# test association of LF2 and NANOG expression
cor.test(cllRNA[thisGeneInd, sharedIDs], latentFactorMat[sharedIDs, "LF2"])
plot(cllRNA[thisGeneInd, sharedIDs], latentFactorMat[sharedIDs, "LF2"])

# also chr12 trisomy and NANOG expression
muts = cllMultiOmics$Mutations
sharedIDs = rnaIDs[rnaIDs %in% colnames(muts)]
trisomyStatus = muts["trisomy12", sharedIDs]
wilcox.test(x = cllRNA[thisGeneInd, sharedIDs][trisomyStatus == 0], y = cllRNA[thisGeneInd, sharedIDs][trisomyStatus == 1])
plot(trisomyStatus, cllRNA[thisGeneInd, sharedIDs])
ggplot(data = data.frame(trisomyStatus, nanogExpr = cllRNA[thisGeneInd, sharedIDs]), 
       mapping = aes(x=trisomyStatus, y=nanogExpr)) + geom_jitter()
#####
# test association of LF2 and NANOG binding region methylation
sharedIDs = colnames(methylMat)[colnames(methylMat) %in% row.names(latentFactorMat)]
meanMethylThisRS=runCOCOA(signal = methylMat, signalCoord = signalCoord, 
                    GRList = GRList[thisRegionSet], 
                    signalCol = colnames(methylMat))
cor.test(latentFactorMat[sharedIDs, "LF2"], as.numeric(meanMethylThisRS[, sharedIDs]))
plot(latentFactorMat[sharedIDs, "LF2"], as.numeric(meanMethylThisRS[, sharedIDs]))
# also test assocation of methylation with chr12 trisomy
sharedIDs = colnames(methylMat)[colnames(methylMat) %in% colnames(muts)]
trisomyStatus = muts["trisomy12", sharedIDs]
wilcox.test(x = as.numeric(meanMethylThisRS[, sharedIDs])[trisomyStatus == 0], y = as.numeric(meanMethylThisRS[, sharedIDs])[trisomyStatus == 1], conf.int = TRUE)

# test whether NANOG expression is associated with methylation in NANOG regions
sharedIDs = rnaIDs[rnaIDs %in% colnames(methylMat)]
plot(as.numeric(meanMethylThisRS[, sharedIDs]), cllRNA[thisGeneInd, sharedIDs]) 
cor.test(as.numeric(meanMethylThisRS[, sharedIDs]), cllRNA[thisGeneInd, sharedIDs], method = "spearman")

# see if MIRA scores are similar
expNanogGRList = GRangesList(resize(GRList[[thisRegionSet]], width=5000, fix="center"))
methylDTList = list()
for (i in 1:ncol(methylMat)) {
    methylDTList[[i]] = as.data.table(cbind(COCOA:::grToDt(signalCoord), methylProp=methylMat[, i]))
}
names(methylDTList) = colnames(methylMat)

miraProfile = lapply(X = methylDTList, FUN = function(x) aggregateMethyl(BSDT = x, 
                                                                         GRList = expNanogGRList, 
                                                                         binNum = 21,
                                                                         minBaseCovPerBin = 0))
nMIRAProfile = rbindNamedList(miraProfile)
miraScores = MIRA::calcMIRAScore(binnedDT = nMIRAProfile)
miraScores = as.data.frame(miraScores)
row.names(miraScores) = miraScores$sampleName
sharedIDs = colnames(methylMat)[colnames(methylMat) %in% row.names(miraScores)]
plot(as.numeric(meanMethylThisRS[, sharedIDs]), miraScores[sharedIDs, ]$score)


sharedIDs = rnaIDs[rnaIDs %in% row.names(miraScores)]
plot(miraScores[sharedIDs, ]$score, cllRNA[thisGeneInd, sharedIDs])
# figure: boxplot of NANOG expression in chr12+,-? correlation plot of NANOG
# expression and LF2?

###############################################################################
# analysis/plots of latent factor 8, WNT signaling

# Wnt3, Wnt5b, Wnt6, Wnt10a, Wnt14, and Wnt16, as well as the Wnt receptor Fzd3, 
# were highly expressed in CLL, compared with normal B cells.
# https://www.ncbi.nlm.nih.gov/pubmed/14973184
myWNTGenes = c(c("WNT3", "WNT5B", "WNT6", "WNT10A", "WNT14", "WNT16"), "WNT5A")
# WNT5A was shown as having a high loading for LF8 in MOFA paper

# ensembl DB
# hg 19
DB = EnsDb.Hsapiens.v75
library(BloodCancerMultiOmics2017)
loadMOFAData(methylMat = TRUE, signalCoord=TRUE, latentFactorMat = TRUE,
             cllMultiOmics=TRUE)

data("dds", package = "BloodCancerMultiOmics2017")
cllRNA = assay(dds)
colnames(cllRNA) <- colData(dds)$PatID
ensemblNames = row.names(cllRNA)
geneDF = genes(filter(DB, filter = GeneNameFilter(value = "WNT", condition = "contains")), return.type="data.frame")
geneDF = dplyr::filter(geneDF, symbol %in% myWNTGenes) 
hist(cllRNA[dplyr::filter(geneDF, symbol == "WNT5A")$gene_id, ])
wnt5aID = dplyr::filter(geneDF, symbol == "WNT5A")$gene_id
MAE = MultiAssayExperiment(list(rna=cllRNA, latentFactors=t(latentFactorMat)))
# shared samples and same ordering
MAE = MAE[, complete.cases(MAE), ]
cor.test(assay(MAE, "rna")[wnt5aID, ], assay(MAE, "latentFactors")["LF8", ])
plot(assay(MAE, "rna")[wnt5aID, ], assay(MAE, "latentFactors")["LF8", ])
# is WNT expression correlated with LF8
mae2 = c(MAE, methyl=methylMat)
mae2 = mae2[, complete.cases(mae2), ]
all(colnames(assay(mae2, "methyl")) == colnames(assay(mae2, "latentFactors")))

# test whether WNT expression is associated with methylation in downstream region sets
thisRegionSet = "E008-H3K4me1.narrowPeak"
thisRegionSet = "GSM1124068_SOX2.bed"
thisRegionSet = "GSM1124071_NANOG.bed"
thisRegionSet = "Pou5f1.bed"
meanMethylThisRS=aggregateSignalGRList(signal = assays(mae2)$methyl, signalCoord = signalCoord, 
                          GRList = GRList[thisRegionSet], 
                          signalCol = colnames(assays(mae2)$methyl))
meanMethylThisRS = meanMethylThisRS[, colnames(assays(mae2)$methyl)]

cor.test(as.numeric(meanMethylThisRS), assays(mae2)$rna[wnt5aID, ])
# for Pou5f1: cor=-0.1484203 , p-value = 0.0858
plot(as.numeric(meanMethylThisRS), assays(mae2)$rna[wnt5aID, ])



cor.test(as.numeric(meanMethylThisRS), assays(mae2)$latentFactors["LF8", ])
plot(as.numeric(meanMethylThisRS), assays(mae2)$latentFactors["LF8", ])
# methylation level is positively correlated with LF8, implying that
# activity is negatively correlated with LF8

cov(as.numeric(meanMethylThisRS), assays(mae2)$latentFactors["LF8", ])
plot(assays(mae2)$latentFactors["LF8", ], assays(mae2)$rna[wnt5aID, ])
cor.test(assays(mae2)$latentFactors["LF8", ], assays(mae2)$rna[wnt5aID, ])
# WNT5a is (significantly) negatively correlated with LF8

# testing other WNT genes
cllMultiOmics$mRNA[geneDF$gene_id, ]
cor.test(cllMultiOmics$mRNA[wnt5aID, colnames(assays(mae2)$latentFactors)], 
                        assays(mae2)$latentFactors["LF8", ])
plot(cllMultiOmics$mRNA[wnt5aID, colnames(assays(mae2)$latentFactors)], 
         sqrt(assays(mae2)$rna[wnt5aID, ]))
cor.test(cllMultiOmics$mRNA[wnt5aID, colnames(assays(mae2)$latentFactors)], 
     as.numeric(meanMethylThisRS))

#############################################################################
# Panel C Meta region profiles
signalCol = lfCols
topRSInd = rsRankingIndex(rsScores = rsScores, signalCol)

topInd = topRSInd[, "LF8"][1:50]
# uTopInd = unique(c(topPC1Ind, topPC2Ind))
uTopInd = topInd

# topRSList = GRList[uTopInd]
topRSList = lapply(X = GRList[uTopInd], FUN = function(x) resize(x = x, width = 14000, fix = "center"))
topRSNames = names(topRSList)



multiProfileP = makeMetaRegionPlots(signal=featureLFCor, 
                                    signalCoord=signalCoord, GRList=topRSList, 
                                    rsNames=topRSNames, 
                                    signalCol=signalCol, binNum=21, aggrMethod ="default") 

ggsave(filename = ffPlot(paste0(plotSubdir,
                                "/metaRegionLoadingProfiles", 
                                inputID, ".pdf")), plot = multiProfileP[["grob"]], device = "pdf", limitsize = FALSE)
# individual mr profiles

# normalize so plots from different PCs will be comparable
multiProfileP2 = list()
multiProfileP2[["metaRegionData"]] = normalizeMRProfile(signal=featureLFCor, signalCol=signalCol, 
                                    multiProfileP$metaRegionData, 
                                    rsNames=names(multiProfileP$metaRegionData), 
                                    normMethod = "zscore")

names(multiProfileP2[["metaRegionData"]])


topRSNames = c("Pou5f1.bed", 
               "GSM1124068_SOX2.bed", 
               "E008-H3K4me1.narrowPeak", "GSM1124071_NANOG.bed")
abbrevNames = c("POU5F1", "SOX2", "H3K4me1_in_H9_cell_line", "NANOG")
# topRSList = topRSList[topRSNames]
# topRSNames = c("GSM835863_EP300.bed", 
#                "GSM607949_GATA1.bed")
minVal = -0.2
maxVal = 0.8
for (i in seq_along(topRSNames)) {
    thisRS = multiProfileP2[["metaRegionData"]][topRSNames[i]]
    # pcP = lapply(X = thisRS, FUN = function(x) tidyr::gather(data = x, key = "PC", value="loading_value", signalCol))
    # pcP = lapply(X = pcP, as.data.table)
    # pcP = lapply(pcP, function(x) x[, PC := factor(PC, levels = signalCol)])
    
    for (j in 8) {
        myPlot = ggplot(data = dplyr::filter(thisRS[[1]], PC %in% signalCol[j]), mapping = aes(x =binID , y = loading_value)) + 
            # ggplot(data = pcP[[1]], mapping = aes(x =binID , y = loading_value)) + 
            geom_line() + ylim(c(minVal, maxVal)) + geom_hline(yintercept = 0, col="red", alpha = 0.25) +
            # facet_wrap(facets = "PC") + 
            ggtitle(label = wrapper(topRSNames[i], width=30)) + xlab("Genome around Region Set, 14 kb") + 
            ylab("Normalized Correlation") + 
            theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
                  axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                  axis.title = element_blank(), title = element_blank()
                  #axis.text.y=element_blank()
            )
        myPlot
        ggsave(filename = ffPlot(paste0(plotSubdir, 
                                        "/mrProfiles", abbrevNames[i],
                                        "_", signalCol[j], ".svg")), 
               plot = myPlot, device = "svg", height=30, width=50, units = "mm")
    }
}

##############################################
# annoScorePlot for LF8
# filter out low-coverage region sets first


.analysisID = dataID
plotUnits = "mm"
plotWidth = 80
plotHeight = 80
stemPattern = paste0(c("POU5F1", "SOX2", "H3K4me1_in_H9_cell_line", "NANOG"), 
                     collapse ="|")

for (i in paste0("LF8")) {
    
    a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = i, 
                          pattern = c(stemPattern), 
                          patternName = c("stem-related")) +
        theme(legend.position = c(0.15, 0.15)) +
        scale_color_manual(values = c("orange", "red")) + 
        xlab(paste0("Region set rank (", i, ")")) + 
        theme(axis.title.y = element_blank(), 
              legend.text = element_blank(), 
              legend.title = element_blank(), 
              legend.position = "none") +
        scale_x_continuous(breaks = c(0, 1000, 2000), 
                           labels= c("0", "1000", "2000"), limits=c(-25, nrow(rsScores) + 25))
    
    a 
    ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", i, "_", .analysisID, ".svg"), 
           plot = a, device = "svg", width = plotWidth, height = plotHeight, units = plotUnits)
}

# making a plot with legend
i="PC1"
a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = i, 
                      pattern = c("esr1|eralpha", "gata3|foxa1|h3r17"), 
                      patternName = c("ER", "ER-related")) +
    theme(legend.position = c(0.15, 0.15)) +
    scale_color_manual(values = c("blue", "red", "orange")) + 
    xlab(paste0("Region set rank (", i, ")"))

a 
ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", i, "_", .analysisID, "_withLegend.svg"), 
       plot = a, device = "svg", width = plotWidth, height = plotHeight, units = plotUnits)


# "filter()" was masked by another package. Get dplyr filter() back
library(dplyr)