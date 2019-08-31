
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(MIRA)

##############################################################################
# Fig. 4, MOFA CLL
# http://msb.embopress.org/content/14/6/e8124
# factors 1 (strong), 7 and 9 were associated with DNA methylation
# factor 1: cell type/differentiation, factor 7: chemo-immunotherapy treatment prior to sample collection
# factor 7: del17p, TP53 mutations, methylation of oncogenes

loadMOFAData(methylMat = TRUE, signalCoord=TRUE, latentFactorMat = TRUE,
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

##############################################################################
# figure connecting LF2 and chr12 trisomy to NANOG
# NANOG is ENSG00000111704
library(BloodCancerMultiOmics2017)
loadMOFAData(methylMat = TRUE, signalCoord=TRUE, latentFactorMat = TRUE,
             cllMultiOmics=TRUE)
data("dds", package = "BloodCancerMultiOmics2017")
cllRNA = assay(dds)
colnames(cllRNA) <- colData(dds)$PatID
head(cllRNA)

thisRegionSet = "GSM1124068_SOX2.bed"
thisGeneInd = grep(pattern = "ENSG00000181449", x = row.names(cllRNA))
# thisRegionSet = "GSM1124071_NANOG.bed"
# thisGeneInd = grep(pattern = "ENSG00000111704", x = row.names(cllRNA))

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
wilcox.test(x = cllRNA[thisGeneInd, sharedIDs], y = trisomyStatus)
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
wilcox.test(x = as.numeric(meanMethylThisRS[, sharedIDs]), y = trisomyStatus, conf.int = TRUE)

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


