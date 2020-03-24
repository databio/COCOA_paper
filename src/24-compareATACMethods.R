# compare COCOA to other methods that work on ATAC-seq data

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
#devtools::load_all(ffCode("COCOA"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "24-compareATACMethods"
plotSubdir = "24-compareATACMethods/"
sheetsDir = ffProc("COCOA_paper/analysis/sheets/")

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 300
removeLowCov=TRUE
covCutoff=100
plotWidth = 100 # mm
plotHeight = 100 # mm
plotUnits = "mm"

######################################################################
# load data

variationMetric = "cor"

# loads signalMat and signalCoord
loadBRCAatac()
genomicSignal = signalMat
simpleCache(paste0("brcaATACPCA_", ncol(genomicSignal)))
# the PC scores
sampleLabels = brcaATACPCA_73$x

dataID = paste0("brcaATAC", ncol(genomicSignal))

colsToAnnotate = paste0("PC", 1:10)

### get shared samples and put data in same order 
sharedSamples = colnames(signalMat)[colnames(signalMat) %in% row.names(sampleLabels)]
genomicSignal = signalMat[, sharedSamples]
sampleLabels = sampleLabels[sharedSamples, colsToAnnotate]

# # loads database of region sets 
# # (assigns GRList, rsName, rsDescription to global environment)
# loadGRList(genomeV="hg38")

# get the "true" COCOA scores before doing the permutation test
simpleCache(paste0("rsScores_", dataID, "_", variationMetric), 
            assignToVariable = "realRSScores")
###########################################################################
# reprocess/reformat data for chromVAR and brockman

# convert peak matrix to RangedSummarizedExperiment for chromVAR
# assay name needs to be "counts" for future function, even though using normalized signal
example_counts = SummarizedExperiment(assays = list(counts=signalMat+abs(min(signalMat))), 
                         rowRanges=signalCoord)#, colData=)
# add coverage info?
#rowData: score, qval, name
# colData: Cell_Type, depth (samples are rows, row.names are the colnames of assay)

#################### for Brockman
# binarize peak matrix to either present or absent
# only necessary for Brockman
# signalMat
######################################################
# @param data. matrix or data.frame?
# @param quantCutoff. Only return data greater than or equal to this quantile
# (fraction)
# @param MARGIN. Filter based on: 1 for rows, 2 for columns. Argument for apply
# @param evalMetric function.
# @param returnLogical logical. If TRUE, return TRUE for indices that are 
# above threshold, FALSE for indices below
# @return data without rows/columns below quantile threshold 
quantileFilter = function(data, quantCutoff=0.8, MARGIN=1, evalMetric=mean, 
                          returnLogical=FALSE) {
    
    cutoffVal = quantile(data, quantCutoff)
    
    if (returnLogical) {
        filtData = data >= cutoffVal
    } else {
        # code not finished yet
        myMetrics = apply(X = data, MARGIN = MARGIN, evalMetric)
        # filtData = data[myMetrics >= cutoffVal]
    }
    
    return(filtData)
}
###########################################################
library(chromVAR)
library(motifmatchr)
library("BSgenome.Hsapiens.UCSC.hg38")

signalMatMin0 = signalMat+abs(min(signalMat))
signalMat50 = quantileFilter(signalMat, quantCutoff=0.5, returnLogical = TRUE)
# could do high confidence peaks, only upper 20% considered open
signalMat80 = quantileFilter(signalMat, quantCutoff=0.8, returnLogical = TRUE)
# peak by sample

# find which motifs are present in each peak
# peak by motif

kmer_ix <- matchKmers(6, example_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
motifs <- getJasparMotifs()
motif_ix <- matchMotifs(motifs, example_counts, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
motif_ix = assays(motif_ix)$motifMatches

# for each sample, count the motifs for the peaks that are present
# result: motif by sample
motifCountsMin0=t(motif_ix) %*% signalMatMin0

# "scaled so that each k-mer had mean 0 and a standard deviation (SD) of 1"

############################################################################
# Chromvar R/Bioconductor package

# BiocManager::install("chromVAR", dependencies=TRUE)
# BiocManager::install("JASPAR2016")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", dependencies=TRUE)
# BiocManager::install("motifmatchr")
library(chromVAR)
library(motifmatchr)
# library(Matrix)
library(SummarizedExperiment)
library("BSgenome.Hsapiens.UCSC.hg38")



# input is aligned fragments

# data(example_counts, package = "chromVAR")
head(example_counts)
head(assay(example_counts, 1))

example_counts <- addGCBias(example_counts, 
                            genome = BSgenome.Hsapiens.UCSC.hg38)
head(rowData(example_counts))

# skip this because data has been preprocessed and we are using uniform peakset
# #find indices of samples to keep
# counts_filtered <- filterSamples(example_counts, min_depth = 1500, 
#                                  min_in_peaks = 0.15, shiny = FALSE)
# counts_filtered <- filterPeaks(counts_filtered, non_overlapping = TRUE)
counts_filtered = example_counts
motifs <- getJasparMotifs()
# out=matches or out=scores can be passed to 
motif_ix <- matchMotifs(motifs, counts_filtered, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
head(assay(motif_ix, 1))
# see matchKmers() to use Kmers instead of motifs
# kmer_ix <- matchKmers(6, counts_filtered, 
#                       genome = BSgenome.Hsapiens.UCSC.hg38)
assayNames(counts_filtered)
dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)

# calculation of final "score" and p-value
variability <- computeVariability(dev)

simpleCache("chromVAR_BRCA_variability", {
    variability
}, assignToVariable = "variability")


################### visualization
# a plot similar to COCOA's region set score distribution
pdf(file = ffPlot(paste0(plotSubdir, "chromVARScoreDist.pdf"))) 
    scorePlot = plotVariability(variability, use_plotly = FALSE) 
    scorePlot
dev.off()
svg(filename = ffPlot(paste0(plotSubdir, "chromVARScoreDist.svg"))) 
    plotVariability(variability, use_plotly = FALSE) 
dev.off()
# do we need to make signal data non negative or transform it somehow? 
# revisit this^

chromScores = variability
colnames(chromScores)[1] = "rsName"

# see 23-atac visualization.R for source, not including GATA3 because it's ER-related
hemaTFs = c("RUNX1", "SCL|TAL1", "PU.1|PU1|SPI1", 
            "CEBPA", "IRF8", "GFI1", "CEBPE", 
            "TCF3", "KLF1", "GATA1", "GATA2", "Ikaros|IKZF1", "CMYB", "NFE2", 
            "TCF3", "EBF1", "PAX5", "FOXO1", "ID2")
hemaTFs = unique(hemaTFs)
hemaPattern = paste0(hemaTFs, collapse = "|")


plotAnnoScoreDist(rsScores = chromScores, colsToPlot = "variability", 
                  pattern = "FOX-family", patternName = "FOX")
plotAnnoScoreDist(rsScores = chromScores, colsToPlot = "variability", 
                  pattern = "ESR|FOXA1|GATA3|H3R17me", patternName = "ER-related")
annoScoreDist = plotAnnoScoreDist(rsScores = chromScores, colsToPlot = "variability", 
                  pattern = c("JUN|FOS|JDP2", "ESR1|FOXA1|GATA3|H3R17me"), 
                  patternName = c("AP1 (JUN or FOS)", "ER-related"),
                  alpha=1) +
theme(legend.position = c(0.6, 0.6)) +
    scale_color_manual(values = c("black", "red", "gray")) + 
    scale_x_continuous(breaks = c(0, 250, 500), 
                       labels= c("0", "250", "500"), limits=c(-25, nrow(chromScores) + 25))

ggsave(filename = ffPlot(paste0(plotSubdir, "chromScoreDistERRelated_withLegend.svg")), 
       plot = annoScoreDist, device = "svg", height = plotHeight, 
       width = plotWidth, units = "mm")
annoScoreDist = annoScoreDist + theme(legend.position = "none")
ggsave(filename = ffPlot(paste0(plotSubdir, "chromScoreDistERRelated_noLegend.svg")), 
       plot = annoScoreDist, device = "svg", height = plotHeight/2, 
       width = plotWidth/2, units = "mm")

annoScoreDist = plotAnnoScoreDist(rsScores = chromScores, colsToPlot = "variability", 
                  pattern = c("FOX", hemaPattern), 
                  patternName = c("FOX family", "Hematopoietic TFs"),
                  alpha=1) +
    theme(legend.position = c(0.6, 0.6)) +
    scale_color_manual(values = c("green", "orange", "gray")) + 
    scale_x_continuous(breaks = c(0, 250, 500), 
                       labels= c("0", "250", "500"), limits=c(-25, nrow(chromScores) + 25))
ggsave(filename = ffPlot(paste0(plotSubdir, "chromScoreDist_FOX_Hema_withLegend.svg")), 
       plot = annoScoreDist, device = "svg", height = plotHeight, 
       width = plotWidth, units = "mm")
annoScoreDist = annoScoreDist + theme(legend.position = "none")
ggsave(filename = ffPlot(paste0(plotSubdir, "chromScoreDist_FOX_Hema_noLegend.svg")), 
       plot = annoScoreDist, device = "svg", height = plotHeight/2, 
       width = plotWidth/2, units = "mm")

View(arrange(variability, desc(variability)))

# point out breast cancer related TFs
# FOXA1 ESR1 ESR2
# JUN/FOS?
# look up other top TFs to see if they are associated with ER


# plot by category
# FOX
# JUN/FOS (JDP2)
# hematopoietic TFs

####################################################################
# BROCKMAN
# https://carldeboer.github.io/brockman.html
# data processing: https://github.com/Carldeboer/Brockman/blob/master/brockman_pipeline
# https://github.com/Carldeboer/AMUSED
# maximum kmer length is 8
# searched sequence for both strands
# k-mer program AMUSED takes fasta as input

# DNA sequences were then extracted from these loci using twoBitToFa [49] and scanned 
# for k-mer content using AMUSED (https://github.com/Carldeboer/AMUSED), considering both DNA strands, to yield a vector of k-mer 
# frequencies for each cell that was used in subsequent analyses, including all gapped k-mers from length 1 to 8. 

library(devtools)
install_github("Carldeboer/BrockmanR")
library("BrockmanR")

# input for BROCKMAN is kmer count matrix (samples x kmers)

# do PCA and tSNE on k-mer matrix
pcs = doKMerPCA(allK562Data, nPCs = "jackstraw");
# here, `pcs` is the object returned by `prcomp`, with several other entries 
# for tSNE and the number of significant PCs


# treatmentPCs = findDistinguishingPCs(pcs$x[,1:pcs$nPCs], pcs$tSNEProj[c("ID","treated")])
# treatmentPCs = treatmentPCs[order(treatmentPCs$P),] 

cellPCProjections = as.data.frame(pcs$x[,1:pcs$nPCs])
cellPCProjections$goodID = row.names(cellPCProjections);
cellPCProjections = merge(cellPCProjections, sampleDesc, by="goodID")


# get enriched TFs for each PC ## if run on all PCs, this can take awhile
#tfEnrichmentsPBM = getKMerTFEnrichment(pcs$rotation[,1:pcs$nPCs], cisbp$binaryPBMZScores); #all PCs
tfEnrichmentsPBM = getKMerTFEnrichment(pcs$rotation[,c(2,4)], cisbp$binaryPBMZScores); # just PCs 2 and 4

# here would be a good place to also consider PWM motifs, and that would look something like the below:
#tfEnrichmentsPWM = getKMerTFEnrichment(pcs$rotation[,c(2,4)], cisbp$binaryPWMScores, n_max = 15000); # just PCs 2 and 4 
# this takes substantially longer because it includes all gapped k-mers, rather than just ungapped kmers (as with cisbp$binaryPBMZScores)
# Because the set size is so much larger, p-values are substantially more significant.  This leads to a lot of false positives because many motifs are so similar to each other.
# We penalized these by adding ln(10^110) to the ln(P) value, effectively making the P-value cutoff of 10^-2 be 10^-112.
# This was based on the high false-positive rate of PWM motifs (which can be highly similar and included a larger number of k-mers),
#  and the elbow of the log(P-value) curves for both PWMs and PBM 8-mer enrichments.

# add TF names to this table
tfEnrichmentsPBM = merge(tfEnrichmentsPBM, cisbp$TFTable[c("Motif_ID","TF_Name")], by="Motif_ID")

# if a TF has both directly-determined motifs and indirectly-determined motifs, consider only those that are directly determined
tfEnrichmentsPBM = preferDirect(tfEnrichmentsPBM, cisbp$TFTable)

#This would also be a good time to filter out TFs that are expressed in your system, but here we will leave all TFs

# Correct for multiple hypothesis testing
tfEnrichmentsPBM$Bon.P = tfEnrichmentsPBM$p + log(nrow(tfEnrichmentsPBM)) +  log(3000) #approximate Bonferroni MHT correction; multiply by 3000 for n_max

# Take significant hits
tfEnrichmentsPBM = tfEnrichmentsPBM[tfEnrichmentsPBM$Bon.P < log(0.01),] # cutoff of P<0.01

# sort by significance
tfEnrichmentsPBM = tfEnrichmentsPBM[order(tfEnrichmentsPBM$p),]

# remove motifs that appear to be redundant, keeping the best motif per PC/enrichment direction, regardless of the corresponding TF
tfEnrichmentsPBM_NR = dropSimilarMotifs(tfEnrichmentsPBM, cisbp$similarMotifs); #similarity defaults to 0.5, which may need to be increased for some applications if too dissimilar motifs are being consolidated

# take best motif per TF, if TFs have more than one motif per PC/enrichment direction
tfEnrichmentsPBM_NR_onePerTF = bestMotifPerTF(tfEnrichmentsPBM_NR);

#sort by significance
tfEnrichmentsPBM_NR_onePerTF = tfEnrichmentsPBM_NR_onePerTF[order(tfEnrichmentsPBM_NR_onePerTF$p),]
head(tfEnrichmentsPBM_NR_onePerTF[tfEnrichmentsPBM_NR_onePerTF$PC!="PC1",1:8]) # PC1 is often highly-correlated with GC content