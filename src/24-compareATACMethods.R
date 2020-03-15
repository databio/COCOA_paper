# compare COCOA to other methods that work on ATAC-seq data

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
devtools::load_all(ffCode("COCOA"))

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

# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg38")

# get the "true" COCOA scores before doing the permutation test
simpleCache(paste0("rsScores_", dataID, "_", variationMetric), 
            assignToVariable = "realRSScores")
###########################################################################
# reprocess data for chromVAR and brockman

# binarize peak matrix to either present or absent
signalMat

# get sequence information for each peak

# for each sample, assign the sample the sequences for peaks that are present

############################################################################
# Chromvar R/Bioconductor package

BiocManager::install("chromVAR", dependencies=TRUE)
BiocManager::install("DirichletMultinomial", dependencies=TRUE)
DirichletMultinomial
library(chromVAR)
library(motifmatchr)
# library(Matrix)
library(SummarizedExperiment)



# input is aligned fragments

data(example_counts, package = "chromVAR")
head(example_counts)

library(BSgenome.Hsapiens.UCSC.hg19)
example_counts <- addGCBias(example_counts, 
                            genome = BSgenome.Hsapiens.UCSC.hg19)
head(rowData(example_counts))

# according to chromvar, each peak must have at least one read for each sample
# therefore, we add a psuedocount

motifs <- getJasparMotifs()
# out=matches or out=scores can be passed to 
motif_ix <- matchMotifs(motifs, counts_filtered, 
                        genome = BSgenome.Hsapiens.UCSC.hg19)
# see matchKmers() to use Kmers instead of motifs
# kmer_ix <- matchKmers(6, counts_filtered, 
#                       genome = BSgenome.Hsapiens.UCSC.hg19)

####################################################################
# BROCKMAN
# https://carldeboer.github.io/brockman.html
# data processing: https://github.com/Carldeboer/Brockman/blob/master/brockman_pipeline
# https://github.com/Carldeboer/AMUSED
# maximum kmer length is 8
# searched sequence for both strands
# k-mer program AMUSED takes fasta as input

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
