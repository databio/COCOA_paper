
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(COCOA)
library("ComplexHeatmap")
library(ggbiplot)
library(readr)

# 
setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)
plotSubdir = "23-brcaATAC/"

if(!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}


# cache data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))
regionCovCutoff = 100
plotWidth = 100 # mm
plotHeight = 100 # mm
plotUnits = "mm"


##########

loadBRCAatac(signalMat = TRUE, signalCoord = TRUE, 
             pcScores = TRUE, loadingMat = TRUE)

# parameters that should have been loaded by previous script
# (load them now if I'm just running this script by itself)
if (!exists("nPerm")) {
    nPerm = 300
}
if (!exists("dataID")) {
    dataID =paste0("brcaATAC", ncol(signalMat))
    # dataID = "brcaATACPCAPerm"
}
if (!exists("variationMetric")) {
    variationMetric = "cor"
}

.analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)
inputID = paste0("_", nPerm, "Perm_", variationMetric, 
                                     "_", "brcaATAC", ncol(signalMat))

simpleCache(paste0("rsScores", paste0("_", "brcaATAC", ncol(signalMat)), 
                   "_", variationMetric), assignToVariable = "rsScores", reload = TRUE)

# simpleCache(paste0("pRankedScores", .analysisID), assignToVariable = "rsScores", reload = TRUE)
# rsEnSortedInd = rsRankingIndex(rsScores = rsScores, 
#                                signalCol = list(paste0(paste0("PC", 1:10), "_PValGroup"), 
#                                                 paste0("PC", 1:10)), 
#                                decreasing = c(TRUE, TRUE), newColName = paste0("PC", 1:10))
simpleCache(paste0("rsPermScores", .analysisID), assignToVariable = "rsPermScores")
signalCol = paste0("PC", 1:10)

loadGRList(genomeV = "hg38")

# screen out region sets that have less than 100 regions that overlap the data
if (length(GRList) != nrow(rsScores)) {
    stop("GRList and rsScores do not match")
}
lowCovRS <- rsScores$regionSetCoverage < 100
rsScores <- rsScores[!lowCovRS, ]
GRList <- GRList[!lowCovRS]
rsName <- rsName[!lowCovRS]
rsDescription <- rsDescription[!lowCovRS]

rsEnSortedInd = rsRankingIndex(rsScores = rsScores, 
                               signalCol = paste0("PC", 1:10), 
                               decreasing = TRUE)

##############################################################################
########## make PCA plot
aMetadata = read.csv(paste0(Sys.getenv("CODE"), "COCOA_paper/metadata/tcga_brca_atacseq_metadata.csv"))
aMetadata$ER_status[aMetadata$ER_status == ""] = NA
pcScores = data.frame(pcScores, pcaNames = row.names(pcScores))
pcScoreAnno= merge(pcScores, aMetadata, by.x = "pcaNames", by.y= "subject_ID", all.x=TRUE)
pcScoreAnno = pcScoreAnno[order(pcScoreAnno$pcaNames)[!duplicated(sort(pcScoreAnno$pcaNames))], ]
# colorClusterPlots(pcScoreAnno, plotCols = paste0("PC", c(1,2)), colorByCols = "ER_status")
pcScoreAnno$ER_status[pcScoreAnno$ER_status == ""] = NA
aPCAPlot = ggplot(data = pcScoreAnno, mapping = aes(x = PC1, y= PC2)) + geom_point(aes(col=ER_status), size = 1, alpha=0.5) +
    coord_fixed() + theme(legend.position = c(0.15, 0.15), axis.text = element_blank(), 
                          axis.ticks = element_blank()) + scale_color_discrete(name="ER status")
aPCAPlot
ggsave(filename = ffPlot(paste0(plotSubdir, "pc1_2_BRCA_ATAC.svg")), 
       height=plotHeight-10, width=plotWidth-10, units = "mm",
       plot = aPCAPlot, device = "svg")
########
# bigger plot for presentations
biggerPCA = aPCAPlot + theme(legend.position = "right",  axis.text = element_text(), axis.ticks = element_line()) + 
    scale_x_continuous(breaks=scales::pretty_breaks(n=3)) + 
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) 
    
biggerPCA

ggsave(filename = ffPlot(paste0(plotSubdir, "pc1_2_BRCA_ATAC_presentation.svg")), 
       height=80, width=100, units = "mm",
       plot = biggerPCA, device = "svg")

# pcaWithAnno = cbind(mPCA$x, patientMetadata[row.names(mPCA$x) ,])

########
# get % variance explained for top PCs
simpleCache("brcaATACPCA_73", assignToVariable = "aPCA")
varExpl = aPCA$sdev^2 / sum(aPCA$sdev^2)
svg(filename = ffPlot(paste0(plotSubdir, "brca_ATAC73_varianceExplained.svg")))
plot(varExpl)
dev.off()



# 0.21017435 0.07013775 0.05563048 0.03980693 0.03129236 0.02875776

#############################################################################
# panel A
# load(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/atac/brca_peak_pca_sample_names.RData")) # pcaNames

hemaTFs = c("RUNX1", "SCL|TAL1", "PU.1|PU1|SPI1", 
            "CEBPA", "IRF8", "GFI1", "CEBPE")

# general: http://www.jbc.org/content/270/10/4955.short
# (table)
# TCF3 (E2A), KLF1 (EKLF), GATA1, GATA2, Ikaros, c-MYB, p45 NF-E2, PAX5, PU.1, RBTN2, SCL/TAL1
hemaTFs = c(hemaTFs, c("TCF3", "KLF1", "GATA1", "GATA2", "Ikaros|IKZF1", "CMYB", "NFE2"))

# lymphoid: https://doi.org/10.1182/blood-2014-12-575688
# (figure 1. table 1) IKZF1, TCF3, EBF1, PAX5, FOXO1, ID2, GATA3
hemaTFs = c(hemaTFs, c("TCF3", "EBF1", "PAX5", "FOXO1", "ID2"))# "GATA3"))
hemaTFs = unique(hemaTFs)
hemaPattern = paste0(hemaTFs, collapse = "|")

myPCs = paste0("PC", 1:4)
for (i in seq_along(myPCs)) {
    pcAnnoScoreDist = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = myPCs[i], 
                                        pattern = c("esr1|eralpha", "foxa1|gata3|H3R17me2", hemaPattern), 
                                        patternName = c("ER", "ER-related", "Hematopoietic TFs")) + 
        theme(legend.position = "none", axis.title.y = element_blank()) +
        scale_color_manual(values = c("blue", "red", "orange", "gray")) + 
        xlab(paste0("Region set rank (", myPCs[i],")")) +
        scale_x_continuous(breaks = c(0, 1000, 2000), 
                           labels= c("0", "1000", "2000"), limits=c(-25, nrow(rsScores) + 25)) + 
        scale_y_continuous(breaks = scales::pretty_breaks(n=3))
    # coord_fixed(ratio = 10)
    pcAnnoScoreDist 
    ggsave(filename = ffPlot(paste0(plotSubdir, myPCs[i], "AnnoScoreDistERRelated.svg")), 
           plot = pcAnnoScoreDist, device = "svg", height = plotHeight/2, width = plotWidth/2, units = "mm")
}

# one with legend
i=1
pcAnnoScoreDist = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = myPCs[i], 
                                    pattern = c("esr1|eralpha", "foxa1|gata3|H3R17me2", hemaPattern), 
                                    patternName = c("ER", "ER-related", "Hematopoietic TFs")) + 
    theme(legend.position = c(0.15, 0.15)) +
    scale_color_manual(values = c("blue", "red", "orange", "gray")) + 
    scale_x_continuous(breaks = c(0, 1000, 2000), 
                       labels= c("0", "1000", "2000"), limits=c(-25, nrow(rsScores) + 25))
# coord_fixed(ratio = 10)
pcAnnoScoreDist 
ggsave(filename = ffPlot(paste0(plotSubdir, myPCs[i], "AnnoScoreDistERRelated_withLegend.svg")), 
       plot = pcAnnoScoreDist, device = "svg", height = plotHeight, width = plotWidth, units = "mm")


# inkscape uses 90 dpi instead of 72
# inkscape total plot size dimensions are 0.8 times ggplot dimensions
pcAnnoScoreDist2 = plotAnnoScoreDist2(rsScores = rsScores, colsToPlot = "PC1", pattern = c("esr|eralpha", "foxa1|gata3|H3R17me2"), 
                   patternName = c("ER", "ER-related"))
pcAnnoScoreDist2 
                                      #pattern = "esr|eralpha|foxa1|gata3|H3R17me2", patternName = "ER-related")
ggsave(filename = ffPlot(paste0(plotSubdir, "pc1AnnoScoreDist2ERRelated.svg")), 
       plot = pcAnnoScoreDist2, device = "svg",  height = 150, width = 150, units = "mm")


################# 
# panel A part 2 
# PC2, immune-related
# reviews of hematopoietic transcription factors (not exhaustive obviously):
# one general review, one myeloid review, and one lymphoid review
# myeloid: https://www.nature.com/articles/nri2024
# (table 1) RUNX1, SCL/TAL1, PU.1, CEBPa, IRF8, GFI1, CEBPe 
# abbrev cebpa
hemaTFs = c("RUNX1", "SCL|TAL1", "PU.1|PU1|SPI1", 
            "CEBPA", "IRF8", "GFI1", "CEBPE")

# general: http://www.jbc.org/content/270/10/4955.short
# (table)
# TCF3 (E2A), KLF1 (EKLF), GATA1, GATA2, Ikaros, c-MYB, p45 NF-E2, PAX5, PU.1, RBTN2, SCL/TAL1
hemaTFs = c(hemaTFs, c("TCF3", "KLF1", "GATA1", "GATA2", "Ikaros|IKZF1", "CMYB", "NFE2"))

# lymphoid: https://doi.org/10.1182/blood-2014-12-575688
# (figure 1. table 1) IKZF1, TCF3, EBF1, PAX5, FOXO1, ID2, GATA3
hemaTFs = c(hemaTFs, c("TCF3", "EBF1", "PAX5", "FOXO1", "ID2")) #, "GATA3"))

hemaTFs = unique(hemaTFs)

hemaPattern = paste0(hemaTFs, collapse = "|")
pc2AnnoScoreDist = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = "PC2", 
                  pattern = hemaPattern, patternName = "Hematopoietic TFs") + 
    theme(legend.position = c(0.15, 0.15)) + scale_color_manual(values = c("red", "orange"))
pc2AnnoScoreDist = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = "PC10", 
                                     pattern = c( "EZH2|SUZ12"),  # hemaPattern,
                                     patternName = c( "polycomb")) + # "Hematopoietic TFs",
    theme(legend.position = c(0.15, 0.15)) + scale_color_manual(values = c("gray", "blue", "red"))

pc2AnnoScoreDist
ggsave(ffPlot(paste0(plotSubdir, "pc2HemaATAC.svg")), 
       plot = pc2AnnoScoreDist, device = "svg", width = plotWidth, 
       height=plotHeight, units = "mm")

pc2AnnoScoreDist = plotAnnoScoreDist2(rsScores = rsScores, colsToPlot = "PC2", 
                   pattern = hemaPattern, patternName = "Hematopoietic TFs")
pc2AnnoScoreDist
ggsave(ffPlot(paste0(plotSubdir, "pc2HemaATAC2.svg")), 
       plot = pc2AnnoScoreDist, device = "svg")

#################
# supplementary fig. immune region sets for other PCs
# supplementary fig. polycomb region sets (EZH2/SUZ12)
if (!dir.exists(ffPlot(plotSubdir, "hemaRSDist/"))) {
    dir.create(ffPlot(plotSubdir, "hemaRSDist/"))
}
pcsToAnnotate = paste0("PC", 1:10)
for (i in seq_along(pcsToAnnotate)) {
    thisPCPlot = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = pcsToAnnotate[i], 
                      pattern = hemaPattern, patternName = "Hematopoietic TFs") + 
        theme(axis.title = element_text(size = 12.5), axis.text = element_text(size=10),
              plot.title = element_text(hjust = 0.5, face="bold")) + 
        scale_color_manual(values = c("red", "lightgray"), guide=FALSE) + 
        ggtitle(label = pcsToAnnotate[i])
    # Title with PC?
    ggsave(ffPlot(paste0(plotSubdir, "hemaRSDist/", pcsToAnnotate[i], "HemaATAC.svg")), 
           plot = thisPCPlot, device = "svg", width = 70, 
           height=70, units = "mm")
    
    # polycomb
    thisPCPlot = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = pcsToAnnotate[i], 
                                   pattern = "EZH2|SUZ12", patternName = "polycomb") + 
        theme(axis.title = element_text(size = 12.5), axis.text = element_text(size=10),
              plot.title = element_text(hjust = 0.5, face="bold")) + 
        scale_color_manual(values = c("lightgray", "red"), guide=FALSE) + 
        ggtitle(label = pcsToAnnotate[i])
    
    ggsave(ffPlot(paste0(plotSubdir, "polycombRSDist_",pcsToAnnotate[i], "_", dataID, ".svg")), 
           plot = thisPCPlot, device = "svg", width = 70, 
           height=70, units = "mm")
    
}

# one with legend

thisPCPlot = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = pcsToAnnotate[i], 
                               pattern = hemaPattern, patternName = "Hematopoietic TFs") + 
    theme(axis.title = element_text(size = 12.5), axis.text = element_text(size=10),
          plot.title = element_text(hjust = 0.5, face="bold")) + 
    scale_color_manual(values = c("red", "lightgray")) + 
    ggtitle(label = pcsToAnnotate[i])
# Title with PC?
ggsave(ffPlot(paste0(plotSubdir, "hemaRSDist/", pcsToAnnotate[i], "HemaATAC_withLegend.svg")), 
       plot = thisPCPlot, device = "svg", width = 70, 
       height=70, units = "mm")

####################
# panel B
# ER status association with PC

# order samples by PC score, color by ER status
erStatusDF = cbind(apply(X = pcScoreAnno[, paste0("PC", 1:4)], 
                         MARGIN = 2, FUN = function(x) order(order(x, decreasing = FALSE))), 
                   as.data.frame(pcScoreAnno)[, "ER_status", drop=FALSE])
erStatusDF = pivot_longer(erStatusDF, cols = c("PC1", "PC2", "PC3", "PC4"), names_to = "PC", values_to = "rank")    

erStatusDF$barHeight = rep(1, nrow(erStatusDF))
erStatusDF = arrange(erStatusDF, PC, rank) 

erStatusPlot = ggplot(data = filter(erStatusDF, PC %in% c("PC1", "PC2")), mapping = aes(x=rank, y=barHeight, group=PC)) + 
    geom_col(aes(col=ER_status)) + scale_color_discrete(breaks=c("Positive", "Negative")) +
    xlab("Samples ordered by PC score") + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                                axis.title.y = element_blank(), axis.line = element_blank())
erStatusPlot


ggplot2::ggsave(filename=ffPlot(paste0(plotSubdir,"/orderedERStatus.svg")), 
                plot = erStatusPlot, device = "svg", height = plotHeight / 2, 
                width = plotWidth, units = plotUnits)
erStatusPlot = ggplot(data = filter(erStatusDF, PC %in% c("PC1", "PC2", "PC3", "PC4")), mapping = aes(x=rank, y=barHeight, group=PC)) + 
    geom_col(aes(col=ER_status)) + scale_color_discrete(breaks=c("Positive", "Negative")) +
    xlab("Samples ordered by PC score") + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                                axis.title.y = element_blank(), axis.line = element_blank())
erStatusPlot

ggplot2::ggsave(filename=ffPlot(paste0(plotSubdir,"/orderedERStatus2.svg")), 
                plot = erStatusPlot, device = "svg", height = plotHeight / 2, 
                width = plotWidth, units = plotUnits)

# correlation of each PC with ER status
pcaWithAnno = as.data.frame(pcScoreAnno)
thesePCs = paste0("PC", 1:4)
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
write.csv(x = resultsDF, file = ffSheets(paste0("pcERStatusRelationship","_", 
                                                dataID, ".csv")),row.names = FALSE)



######### make "meta-region" loading profiles
# panel C
atacCor = cor(x = t(signalMat), y = pcScores[, signalCol])
all(colnames(signalMat) == row.names(pcScores))

# load region sets
# source(ffProjCode("load_process_regions_brca.R"))

topPC1Ind = rsEnSortedInd[, "PC1"][1:15]
topPC2Ind = rsEnSortedInd[, "PC2"][1:15]
uTopInd = unique(c(topPC1Ind, topPC2Ind))

# topRSList = GRList[uTopInd]
topRSNames = as.character(rsScores$rsName[uTopInd])

topRSNames = unique(c(topRSNames, 
                      "wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak", 
                      "Human_MCF-7_ESR1_E2-6hr_Jin.bed", 
                      "GSM1501162_CEBPA.bed", "GSM1097879_ERG.bed",
                      "Human_MCF-7_H3R17me2_No-treatment_Brown.bed", 
                      "Human_MCF-7_FoxA1_No-treatment_Brown.bed",
                      "wgEncodeAwgTfbsSydhH1hescSuz12UcdUniPk.narrowPeak",
                      "E104-H3K27me3.narrowPeak",
                      "E032-H3K9me3.narrowPeak", 
                      "wgEncodeAwgTfbsBroadH1hescEzh239875UniPk.narrowPeak",
                      "wgEncodeAwgTfbsBroadGm12878Ezh239875UniPk.narrowPeak"))
# "wgEncodeAwgTfbsBroadGm12878Ezh239875UniPk.narrowPeak" top polycomb region set for PC4

topRSList = lapply(X = GRList[topRSNames], FUN = function(x) resize(x = x, width = 14000, fix = "center"))






multiProfileP = makeMetaRegionPlots(signal=atacCor, 
                                    signalCoord=signalCoord, GRList=topRSList, 
                                    rsNames=topRSNames, 
                                    signalCol=signalCol, binNum=21, aggrMethod ="simpleMean") 
multiProfileP2 = makeMetaRegionPlots(signal=atacCor, 
                                    signalCoord=signalCoord, GRList=topRSList, 
                                    rsNames=topRSNames, 
                                    signalCol=signalCol, binNum=21, aggrMethod ="proportionWeightedMean", 
                                    absVal = TRUE, normMethod = "none") 

ggsave(filename = ffPlot(paste0(plotSubdir,
                         "/metaRegionLoadingProfiles", 
                         inputID, ".pdf")), plot = multiProfileP[["grob"]], device = "pdf", limitsize = FALSE)
ggsave(filename = ffPlot(paste0(plotSubdir, "/metaRegionLoadingProfilesWeightedMean", inputID, ".pdf")), 
       plot = multiProfileP2[["grob"]], device = "pdf", limitsize = FALSE)
# individual mr profiles

# normalize
multiProfileP2 = normalizeMRProfile(signal=atacCor, signalCol=signalCol, 
                                    multiProfileP2$metaRegionData, 
                                    names(multiProfileP2$metaRegionData),
                                    normMethod = "zscore")

multiProfileP2[["metaRegionData"]]=multiProfileP2
names(multiProfileP2[["metaRegionData"]])

topRSNames = c("wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak", 
               "Human_MCF-7_ESR1_E2-6hr_Jin.bed", 
               "GSM1501162_CEBPA.bed", "GSM1097879_ERG.bed",
               "Human_MCF-7_H3R17me2_No-treatment_Brown.bed", 
               "Human_MCF-7_FoxA1_No-treatment_Brown.bed",
               "wgEncodeAwgTfbsSydhH1hescSuz12UcdUniPk.narrowPeak",
               "E104-H3K27me3.narrowPeak",
               "E032-H3K9me3.narrowPeak", 
               "wgEncodeAwgTfbsBroadH1hescEzh239875UniPk.narrowPeak",
               "wgEncodeAwgTfbsBroadGm12878Ezh239875UniPk.narrowPeak")
abbrevNames = c("GATA3", "ESR1", "CEBPA", "ERG", 
                "H3R17me2", "FOXA1",
                "SUZ12", "H3K27me3", "H3K9me3",
                "EZH2", "EZH2_topRS_PC4")

# topRSNames = c("GSM835863_EP300.bed", 
#                "GSM607949_GATA1.bed")
minVal = -1
maxVal = 2
pcP = list()
for (i in seq_along(topRSNames)) {
    thisRS = multiProfileP2[["metaRegionData"]][topRSNames[i]]
    pcP[1] = thisRS
    # pcP = lapply(X = thisRS, FUN = function(x) tidyr::gather(data = x, key = "PC", value="loading_value", signalCol))
    # pcP = lapply(X = pcP, as.data.table)
    # pcP = lapply(pcP, function(x) x[, PC := factor(PC, levels = signalCol)])
    
    for (j in seq_along(signalCol[1:4])) {
        myPlot = ggplot(data = filter(pcP[[1]], PC %in% signalCol[j]), mapping = aes(x =binID , y = loading_value)) + 
            # ggplot(data = pcP[[1]], mapping = aes(x =binID , y = loading_value)) + 
            geom_line() + ylim(c(minVal, maxVal)) + geom_hline(yintercept = 0, col="red", alpha = 0.25) +
            # facet_wrap(facets = "PC") + 
            ggtitle(label = wrapper(topRSNames[i], width=30)) + xlab("Genome around Region Set, 14 kb") + 
            ylab("Normalized Correlation") + 
            theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
                  axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                  axis.title = element_blank(), title = element_blank(), 
                  axis.text.y=element_blank()
                  )
        myPlot
        ggsave(filename = ffPlot(paste0(plotSubdir, 
                                        "/mrProfilesWeightedMean_", abbrevNames[i],
                                        "_", signalCol[j], ".svg")), 
               plot = myPlot, device = "svg", height=20, width=30, units = "mm")
    }

    
}


##############################################################################
# supplementary figure: raw data in ER regions, ordered by PC score

topRSInd = rsRankingIndex(rsScores = rsScores,
                               signalCol = paste0("PC", 1:10), 
                               decreasing = TRUE, 
                          newColName = paste0("PC", 1:10))
topRSNames = names(GRList)[unique(unlist(topRSInd[1:5, paste0("PC", 1:5)]))]
# list(paste0(paste0("PC", 1:10), "_PValGroup")
abbrevNames = gsub(pattern = ".", replacement = "_", x = topRSNames, fixed = TRUE)

# top 2 region sets from PC1 and PC2
# wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak
# Human_MCF-7_ESR1_E2-6hr_Jin.bed
# GSM1501162_CEBPA.bed
# GSM1097879_ERG.bed

# get uniform scale 
max(signalMat)
hist(signalMat)
apply(X = signalMat, MARGIN = 2, FUN = median)

# scaleMax = max(signalMat)
# scaleMin = min(signalMat)
# aMean = mean(signalMat)
# aSD = sd(signalMat)
# threeSDUp = aMean + 3 * aSD
# threeSDDown = aMean - 3 * aSD
# scaleMax = min(threeSDUp, scaleMax)
# scaleMin = max(threeSDDown, scaleMin)
# scaleMid = (scaleMin + scaleMax) / 2
quantileMat = matrix(ecdf(signalMat)(signalMat), nrow = nrow(signalMat), ncol=ncol(signalMat))
scaleMin = 0
scaleMid = 0.5
scaleMax = 1

myPCs = paste0("PC", 1:4)
colnames(quantileMat) = gsub(pattern = "-",
                                        replacement = "_", 
                                        x = colnames(signalMat))
row.names(aPCA$x) = gsub(pattern = "-",
                         replacement = "_", 
                         x = row.names(aPCA$x))
for (i in seq_along(topRSNames)) {
    thisRSM = COCOA:::averagePerRegion(signal=quantileMat, 
                                       signalCoord = signalCoord, 
                                       regionSet = GRList[[topRSNames[i]]], 
                                       signalCol = colnames(quantileMat))
    thisRSCovScores = COCOA:::averagePerRegion(signal=loadingMat[, myPCs], 
                                               signalCoord = signalCoord, 
                                               regionSet = GRList[[topRSNames[i]]], 
                                               signalCol = myPCs)
    thisRSM = as.data.frame(thisRSM)
    thisRSCovScores = as.data.frame(thisRSCovScores)
    for (j in seq_along(myPCs)) {
        
        # pdf(file = ffPlot(paste0(plotSubdir, 
        #                          "methylAlong", myPCs[j], "_", abbrevNames[i], ".pdf")))
        png(filename = ffPlot(paste0(plotSubdir,
                                     "signalAlong", myPCs[j], "_", abbrevNames[i], ".png")),
            width = 100, height = 75, units = "mm", pointsize = 12,
            bg = "white",  res = 300,
            type = c("cairo", "cairo-png", "Xlib", "quartz"))
        # # svg(filename = ffPlot(paste0(plotSubdir, "methylAlong", myPCs[j], "_", abbrevNames[i], ".svg")))
        draw(signalAlongAxis(genomicSignal = thisRSM[, !(colnames(thisRSM) %in% c("chr", "start", "end"))], 
                             signalCoord = thisRSM[, c("chr", "start", "end")], 
                             regionSet = GRList[[topRSNames[i]]], 
                             sampleScores= aPCA$x[, myPCs], orderByCol=myPCs[j], 
                             topXVariables=100, variableScores = thisRSCovScores[, myPCs[j]], 
                             cluster_columns=TRUE, show_row_names=FALSE, show_column_names=FALSE, 
                             column_title = "Regions (top 100)", name = "chromatin accessibility (read counts)", 
                             row_title = "Samples ordered by PC score", 
                             row_title_gp = gpar(fontsize = 14), # 54
                             column_title_gp = gpar(fontsize = 14),
                             show_heatmap_legend = FALSE, 
                             col=circlize::colorRamp2(breaks = c(scaleMin,scaleMid,scaleMax), 
                                                      colors = c("blue", "#EEEEEE", "red"))
        )
        ) 
        dev.off()
        
        if ((j == 1) && (i == 1)) {
            svg(filename = ffPlot(paste0(plotSubdir, "signalAlong", myPCs[j], "_", abbrevNames[i], "_hasLegend.svg")), 
                width = 10, height=10) # inches
            draw(signalAlongAxis(genomicSignal = thisRSM[, !(colnames(thisRSM) %in% c("chr", "start", "end"))], 
                                 signalCoord = thisRSM[, c("chr", "start", "end")], 
                                 regionSet = GRList[[topRSNames[i]]], 
                                 sampleScores= aPCA$x[, myPCs], orderByCol=myPCs[j], 
                                 topXVariables=10, variableScores = thisRSCovScores[, myPCs[j]], 
                                 cluster_columns=TRUE, show_row_names=FALSE, show_column_names=FALSE, 
                                 column_title = "Regions (top 100)", name = "chromatin accessibility (read counts)", 
                                 row_title = "Samples ordered by PC score",
                                 col=circlize::colorRamp2(breaks = c(scaleMin,scaleMid,scaleMax), 
                                                          colors = c("blue", "#EEEEEE", "red")))) 
            dev.off()
            # row_title_gp = gpar(fontsize = 14), # 54
            # column_title_gp = gpar(fontsize = 14))
        }
        
    }
}



##########################################################################
# does PC2 correspond to immune signature?
mae = curatedTCGAData(diseaseCode = "BRCA", "mutation", dry.run = FALSE)
immuneInfo = colData(mae)[, c("patient.samples.sample.portions.portion.slides.slide.percent_lymphocyte_infiltration",      
"patient.samples.sample.portions.portion.slides.slide.percent_monocyte_infiltration",        
"patient.samples.sample.portions.portion.slides.slide.percent_neutrophil_infiltration")] 
names(immuneInfo) = c("lymphocyte", "monocyte", "neutrophil")
immuneInfo = immuneInfo[row.names(pcScores), ]
dim(immuneInfo)
dim(pcScores)
cor.test(immuneInfo$lymphocyte, pcScores$PC2)
cor.test(immuneInfo$monocyte, pcScores$PC2)
cor.test(immuneInfo$neutrophil, pcScores$PC2)
plot(immuneInfo$monocyte, pcScores$PC2)
plot(immuneInfo$lymphocyte, pcScores$PC2)

############################################################################
rsEnrichment = rsScores
coordinateDT = COCOA::grToDt(signalCoord)
mPCA = list()
mPCA$rotation = loadingMat
mPCA$x = pcScores
mPCA$center = rowMeans(methylData)
methylData = signalMat

PCSTOANNOTATE = paste0("PC", 1:10)


### plots that will be created and script specific parameters for them 
# "comparePCHeatmap"
PCsToAnnotate_cPCH = PCSTOANNOTATE
# "methylAlongPC"
topRSToPlotNum = 15
PCsToAnnotate_mAPC = PCSTOANNOTATE[1:10]
# "regionQuantileByPC"
PCsToAnnotate_rQBPC = PCSTOANNOTATE
topRSInd_rQBPC = unique(unlist(rsEnSortedInd[1:15, ])) # get top region sets from each PC
# pcFromSubset Correlation Heatmap
PCsToAnnotate_pcFSCH = PCSTOANNOTATE
topRSInd_pcFSCH = unique(unlist(rsEnSortedInd[1:15, ])) # get top region sets from each 
## "region set Overlapping Cytosine Proportion" (rsOLCP)
## proportion of cytosines from region set that are shared with other region set
topRSInd_rsOLCP = unique(unlist(rsEnSortedInd[1:10, ]))
## "meta region loading profiles" (mrLP)
topRSInd_mrLP = unique(unlist(rsEnSortedInd[1:10, ]))
PCsToAnnotate_mrLP = PCSTOANNOTATE

# the pipeline
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/COCOA_vis_pipeline.R")) 

