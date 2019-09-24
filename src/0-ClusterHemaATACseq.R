# process hematopoietic cell ATACseq data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74912
# hg19
# to view on UCSC, first go to: 
# https://genome.ucsc.edu/cgi-bin/hgHubConnect?hgsid=752998689_zIkSMeQ40MJXa4la5Iu3gOMjZdXo
# then paste the below URL into the box on the My Hubs page
# https://s3-us-west-1.amazonaws.com/chang-public-data/2016_NatGen_ATAC-AML/hub.txt


library(data.table)
library(cluster)
library(ggplot2)
library(tidyr) # gather, filter?
library(GenomicRanges)
library(factoextra) # for analysis of cluster number
library(simpleCache)
library(ComplexHeatmap)
library(projectInit)
library("BSgenome.Hsapiens.UCSC.hg19") # get sequence, calculate GC content
library(biovizBase) # GCcontent
library(cqn) # GC normalization
library(preprocessCore) # quantile normalization

# runs 00-init.R and sets up environment
projectInit(codeDir = paste0(Sys.getenv("CODE"), "aml_e3999/"), procDir = paste0(Sys.getenv("PROCESSED"), "/aml_e3999/")) 
source(paste0(Sys.getenv("CODE"), "aml_e3999/src/00-init.R"))

Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "aml_e3999/analysis/plots/"))
plotSubdir = "0-hemaAtac/"
Sys.setenv("PLOT.SUBDIR"=plotSubdir)
scriptID = "0hemaAtac"
# saving plots in script subdirectory in plots directory, if not present then create it now
if (!dir.exists(paste0(Sys.getenv("PLOTS"), plotSubdir))) {
    dir.create(paste0(Sys.getenv("PLOTS"), plotSubdir), recursive = TRUE)
}


setwd(paste0(Sys.getenv("PROCESSED"), "/aml_e3999/prjResources/"))
# setwd("C:/Users/John L/Dropbox/Research_Files/data/GEOdata")


######
regionCounts = fread("GSE74912_ATACseq_All_Counts.txt")
regionString = paste0(regionCounts$Chr, "_", regionCounts$Start, "_", regionCounts$End)
regionCounts = regionCounts[, 4:ncol(regionCounts)]
sample_names = colnames(regionCounts)

# to simplify visualization
# regionCounts[regionCounts > 100] = 100
hist(regionCounts$`4983-1A`)
boxplot(regionCounts$`4983-1A`)
sampleMean = colMeans(regionCounts)
hist(sampleMean)
sampleMedian = apply(regionCounts, 2, median)
hist(sampleMedian)

# function for making GenomicRanges object from regionString
# (string with genomic coordinates)
regionStringToGR <- function(regionString) {
    if (length(regionString) == 0) {
        return(GRanges())
    }
    
    chr = sub(pattern = "_.+_.+", replacement = "", regionString)
    start = sub(pattern = "chr.{1,2}_", "", regionString)
    start = as.numeric(sub(pattern = "_.+$", "", start))
    end = as.numeric(sub(pattern = ".+_.+_", "", regionString))
    regionGR = MIRA:::dtToGr(data.frame(chr, start, end))
    return(regionGR)
}

allRegionGR = regionStringToGR(regionString)

#############################################################################
# normalize according to description in original paper
# quantile normalize
regionCountsUnNorm = regionCounts
regionCounts = as.data.frame(normalize.quantiles(as.matrix(regionCounts))) # treats each column as a sample
colnames(regionCounts) = colnames(regionCountsUnNorm)

# get GC content for regions and normalize with cqn package
hg19Genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
gcPerRegion = biovizBase::GCcontent(hg19Genome, allRegionGR)
nRegions = nrow(regionCounts)
nSamples = ncol(regionCounts)
cqnNorm = cqn(as.data.frame(regionCounts)[1:nRegions, 1:nSamples],
                  lengths = rep(501, nRegions), x = gcPerRegion[1:nRegions],
                sizeFactors = colSums(regionCounts)[1:nSamples], 
                verbose = FALSE, lengthMethod = "fixed")
# from docs: Final corrected values are equal to value$y + value$offset.
normCounts = as.data.frame(cqnNorm$y + cqnNorm$offset)

# cqn.subset <- cqn(as.data.frame(regionCounts)[1:nrow(regionCounts), 1:132], lengths = rep(501, 590650), 
#                                    x = gcPerRegion[1:590650], sizeFactors = colSums(regionCounts)[1:132],
#                                     verbose = FALSE, lengthMethod = "fixed")
# should regions be 500 or 501?

##############################################################################
# getting regions accessible in any sample (general hematopoietic open chromatin regions, 
# include cancer?)

# selecting top 15% of regions (highest scores), ~88600 regions
# using 85th percentile as cutoff so some may have more regions if other percentiles
# share the same number (percentile .85 may be same as percentile .84)

# selecting 
sampleQuant = apply(normCounts, 2, function(x) quantile(x, 0.85))
# creating filtering mask based on 85th p
screenMat = mapply(function(x, y) x >= y, normCounts, sampleQuant)
presentInAny = apply(X = screenMat, MARGIN = 1, FUN = any)

# for each sample, selecting regions that have highest scores
sampleRegionGRList = list()
for (i in seq_along(sample_names)) {
    sampleRegionGRList[[i]] = regionStringToGR(regionString[screenMat[, i]])
}
names(sampleRegionGRList) = names(normCounts)
sampleRegionGRList = GRangesList(sampleRegionGRList)

# GRangesList with a separate region set for each sample
save(sampleRegionGRList, file="GSE74912_Region_Sets.RData")

# GRanges object with regions that are accessible in any sample
bloodAccessibleRegions = regionStringToGR(regionString[presentInAny])
save(bloodAccessibleRegions, file="bloodAccessibleRegions.RData")



############# creating consensus region sets #############################
# creating a single region set for each cell type by
# merging count data by cell type then narrowing down to top regions

# metadata created by geofetch.py
metaDT = fread("annocomb_GSE74912.csv")

# fixing some inconsistencies in metadata file
metaDT$Sample_title = sub("\\.trim", replacement = "", metaDT$Sample_title) # 3 to fix
metaDT$Sample_title = sub(" ", "-", metaDT$Sample_title) # "SU575 pHSC" is only target
# adding a sample that was missing from metadata but included in processed count matrix
# all Blast samples are Leukemia and vice versa so I think it is safe to assume that these are Blasts
blankRow = rep("", ncol(metaDT))
names(blankRow) = colnames(metaDT)
# convert back to just data.frame because data.table was giving problems
metaDT = rbind(as.data.frame(metaDT, stringsAsFactors=FALSE), blankRow)
metaDT[nrow(metaDT), "Sample_title"] = "SU072-Leuk-Rep1-150212"
metaDT[nrow(metaDT), "Sample_source_name_ch1"] = "Blast"
metaDT = rbind(metaDT, blankRow)
metaDT[nrow(metaDT), "Sample_title"] = "SU072-Leuk-Rep2-150212"
metaDT[nrow(metaDT), "Sample_source_name_ch1"] = "Blast"

## this code was used before I quantile normalized and GC corrected.
# # first normalize each sample for total reads so samples will be more comparable
# # assumes that each sample that will be combined 
# # (samples from same patient and then later samples from same cell types)
# # has a similar amount of open chromatin
# sampleCount = colSums(regionCounts)
# max(sampleCount)
# min(sampleCount)
# hist(sampleCount) # over an order of magnitude difference (almost 2)
# # arbitrary factor of 600,000 so not dealing with as many fractions (min(sampleCount)=643360)
# # (percentile will be used later so exact numbers don't matter)
# normCounts = mapply(FUN = function(x, y) as.matrix(regionCounts)[, x] * (600000 / y), x=1:ncol(regionCounts), y=sampleCount)
# normCounts = as.data.table(normCounts)
# colnames(normCounts) = colnames(regionCounts)


# merging counts from the same patient so all patients will be weighted the same
donorid = gsub(pattern = "-.+", replacement = "", metaDT$Sample_title)
donorCell = paste0(donorid, "-", metaDT$Sample_source_name_ch1)
uDonorCell = unique(donorCell)
metaDT = cbind(metaDT, donorCell)
mergedCounts = data.frame(matrix(NA, nrow=nrow(normCounts), ncol=length(uDonorCell)))
rownames(mergedCounts) = regionString
colnames(mergedCounts) = uDonorCell
# combine counts if multiple for a single donor/cell type combination
for (i in seq_along(uDonorCell)) {
    # get name of all samples that should be merged
    sharedDonorCell = metaDT[donorCell == uDonorCell[i], "Sample_title"] 
    mergedCounts[, i] = rowMeans(normCounts[, sharedDonorCell, drop=FALSE])
}

cellType = unique(metaDT$Sample_source_name_ch1) # 18 cell types

####
# making new count matrix, merged by cell type (mean among patients)
cellTypeCounts = matrix(nrow = nrow(normCounts), ncol = length(cellType))
# ordering metadata table according to table of counts
# for metadata samples are rows, for count matrix samples are columns
### row.names(metaDT) = metaDT$Sample_title
### ordMetaDT = metaDT[order(colnames(regionCounts)), ] 
# get consensus (sum) of counts for each cell type
for (i in seq_along(cellType)) {
    cellInd = grep(pattern = cellType[i], colnames(mergedCounts))
    # cellInd = ordMetaDT$Sample_source_name_ch1 == cellType[i]
    # average all counts for this cell type by region
    if (length(cellInd) > 1) {
        consensus = rowMeans(mergedCounts[, cellInd])
    } else {
        # "CD34 Cord Blood" and "CD34 Bone Marrow" 
        consensus = mergedCounts[, cellInd]
    }
        
        cellTypeCounts[, i] = consensus
}
colnames(cellTypeCounts) = cellType

###############
# once again taking the top regions 

sampleQuant = apply(cellTypeCounts, 2, function(x) quantile(x, 0.85))
# creating filtering mask based on Xth percentile
screenMat = mapply(function(x, y) x >= y, as.data.frame(cellTypeCounts), sampleQuant)

# for each sample, selecting regions that have highest scores
ctRegionGRList = list()
for (i in seq_along(cellType)) {
    ctRegionGRList[[i]] = regionStringToGR(regionString[screenMat[, i]])
}
names(ctRegionGRList) = cellType
ctRegionGRList = GRangesList(ctRegionGRList)

# GRangesList with a separate region set for each cell type
# top regions but only in respect to the cell type itself, not to other cell types
save(ctRegionGRList, file="GSE74912_Cell_Type_Consensus_Region_Sets.RData")

# select regions that are associated with a single cell type
# keep regions that are high (open) in one cell type but not in others
# start with unfiltered data (has been merged by patient and cell type)


##############################################################################
# Get ATAC region sets in a supervised way
# see figure1 at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5042844/ for
# hematopoietic tree
# 
# after getting each region set, merge close region together (gap of 200 or less)
mergeGapLength = 200 + 1 # any gaps less than this length are merged so add 1

# see code section above for why I'm taking out certain cell types (not in hematopoietic tree)
mainCellMat = cellTypeCounts[, !(colnames(cellTypeCounts) %in% c("CD34 Cord Blood", 
                                                                "CD34 Bone Marrow", "LSC", "Blast", "pHSC"))]
# arbitrarily calling regions above the Xth percentile as being open
cellQuant = apply(mainCellMat, 2, function(x) quantile(x, 0.9))
# calling regions below Xth percentile low/closed
cellQuantLow = apply(mainCellMat, 2, function(x) quantile(x, 0.5))
# creating filtering mask based on Xth percentile
openCellMat = mapply(function(x, y) x >= y, as.data.frame(mainCellMat), cellQuant)
lowCellMat = mapply(function(x, y) x < y, as.data.frame(mainCellMat), cellQuantLow)

# function for getting specific overlaps
# @param presentMat, rows are regions, columns are cell types. TRUE if region present,
# FALSE if absent
# @param presentOR present in any of these cell types eg cell1 or cell2
# @param notPresentOR not present in any of these cell types eg !(cell1 | cell2)
# @param notPresentMat matrix/data.frame. rows are regions, columns are cell types
# same ordering as presentMat. TRUE if region is not present/ is low in these cell types
# @value a boolean with TRUE values for desired regions (rows of presentMat)
subtypeRegionSet = function(presentMat, presentOR, notPresentOR, notPresentMat=NULL) {
    if (!all(presentOR %in% colnames(presentMat))) {
        stop("check presentOR names.")
    }
    if (!all(notPresentOR %in% colnames(presentMat))) {
        stop("check notPresentOR names.")
    }
    
    # desired regions
    presOR = apply(X = as.matrix(presentMat[, presentOR]), MARGIN = 1, any)
    if (!is.null(notPresentMat)) {
        # not present in all these cell types
        nPresOR = apply(X = as.matrix(notPresentMat[, notPresentOR]), MARGIN = 1, all) 
    } else {
        # not present in any of these cell types
        nPresOR = !apply(X = as.matrix(presentMat[, notPresentOR]), MARGIN = 1, any)
    }
    
    # present in desired cell types but not in undesired cell types
    rSubset = presOR & nPresOR
    return(rSubset)
}
# test = matrix(c(T, T, T, F, F, T, F, T), nrow=4)
# sum(test[, 1] & !test[, 2]) 
# sum(!test[, 1] & test[, 2]) 
# colnames(openCellMat)
# "HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "Mono", "Bcell", "CD4Tcell", "CD8Tcell", "NKcell", "CLP", "Ery"     
# missing granulocyte and megakaryocyte
# Present in HSC but not terminally differentiated cells
# Present in terminally differentiated cells but not HSCs
hscNotDiff = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                              presentOR = c("HSC"), 
                              notPresentOR = c("Bcell", "CD4Tcell", "CD8Tcell", "NKcell", "Mono", "Ery"))]
hscNotDiffGR = regionStringToGR(hscNotDiff)
save(hscNotDiffGR, file="hscNotDiffNotMergedGR.RData")
hscNotDiffMergedGR = reduce(hscNotDiffGR, min.gapwidth=mergeGapLength)
save(hscNotDiffMergedGR, file="hscNotDiffGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(hscNotDiffMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/hscNotDiff.bed")

diffNotHSC = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                           presentOR = c("Bcell", "CD4Tcell", "CD8Tcell", "NKcell", "Mono", "Ery"), 
                                           notPresentOR = c("HSC"))]
diffNotHSCGR = regionStringToGR(diffNotHSC)
save(diffNotHSCGR, file="diffNotHSCNotMergedGR.RData")
diffNotHSCMergedGR = reduce(diffNotHSCGR, min.gapwidth=mergeGapLength)
save(diffNotHSCMergedGR, file="diffNotHSCGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(diffNotHSCMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/diffNotHSC.bed")

# Present in HSC or MPP but not terminally differentiated cells
stemNotDiff = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                            presentOR = c("HSC", "MPP"), 
                                            notPresentOR = c("Bcell", "CD4Tcell", "CD8Tcell", "NKcell", "Mono", "Ery"))]
stemNotDiffGR = regionStringToGR(stemNotDiff)
save(stemNotDiffGR, file="stemNotDiffNotMergedGR.RData")
stemNotDiffMergedGR = reduce(stemNotDiffGR, min.gapwidth=mergeGapLength)
save(stemNotDiffMergedGR, file="stemNotDiffGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(stemNotDiffMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/stemNotDiff.bed")


diffNotStem = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                           presentOR = c("Bcell", "CD4Tcell", "CD8Tcell", "NKcell", "Mono", "Ery"), 
                                           notPresentOR = c("HSC", "MPP"))]
diffNotStemGR = regionStringToGR(diffNotStem)
save(diffNotStemGR, file="diffNotStemNotMergedGR.RData")
diffNotStemMergedGR = reduce(diffNotStemGR, min.gapwidth=mergeGapLength)
save(diffNotStemMergedGR, file="diffNotStemGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(diffNotStemMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/diffNotStem.bed")




# present in any myeloid terminally diff. cells but not in any lymphoid diff. cells
# lymphoid but not myeloid
myeloNotLymph = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                              presentOR = c("Mono", "Ery"), 
                                              notPresentOR = c("Bcell", "CD4Tcell", "CD8Tcell", "NKcell"))]
myeloNotLymphGR = regionStringToGR(myeloNotLymph)
save(myeloNotLymphGR, file="myeloNotLymphNotMergedGR.RData")
myeloNotLymphMergedGR = reduce(myeloNotLymphGR, min.gapwidth=mergeGapLength)
save(myeloNotLymphMergedGR, file="myeloNotLymphGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(myeloNotLymphMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/myeloNotLymph.bed")

lymphNotMyelo = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                              presentOR = c("Bcell", "CD4Tcell", "CD8Tcell", "NKcell"), 
                                              notPresentOR = c("Mono", "Ery"))]
lymphNotMyeloGR = regionStringToGR(lymphNotMyelo)
save(lymphNotMyeloGR, file="lymphNotMyeloNotMergedGR.RData")
lymphNotMyeloMergedGR = reduce(lymphNotMyeloGR, min.gapwidth=mergeGapLength)
save(lymphNotMyeloMergedGR, file="lymphNotMyeloGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(lymphNotMyeloMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/lymphNotMyelo.bed")

# compare myeloid branches
# present in monocyte (gran. not present in dataset) but not in erythrocyte (megakar. not present in dataset)
# also erythr. but not monocyte
monoNotEryRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, presentOR = "Mono", notPresentOR = "Ery")]
monoNotEryGR = regionStringToGR(monoNotEryRegions)
save(monoNotEryGR, file="monoNotEryNotMergedGR.RData")
monoNotEryMergedGR = reduce(monoNotEryGR, min.gapwidth=mergeGapLength)
save(monoNotEryMergedGR, file="monoNotEryGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(monoNotEryMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/monoNotEry.bed")

eryNotMonoRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, presentOR = "Ery", notPresentOR = "Mono")]
eryNotMonoGR = regionStringToGR(eryNotMonoRegions)
save(eryNotMonoGR, file="eryNotMonoNotMergedGR.RData")
eryNotMonoMergedGR = reduce(eryNotMonoGR, min.gapwidth=mergeGapLength)
save(eryNotMonoMergedGR, file="eryNotMonoGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(eryNotMonoMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/eryNotMono.bed")

# T cells vs B cells and vice versa
tCellNotBCellRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                                     presentOR = c("CD4Tcell", "CD8Tcell"), 
                                                     notPresentOR = "Bcell")]
tCellNotBCellGR = regionStringToGR(tCellNotBCellRegions)
save(tCellNotBCellGR, file="tCellNotBCellNotMergedGR.RData")
tCellNotBCellMergedGR = reduce(tCellNotBCellGR, min.gapwidth=mergeGapLength)
save(tCellNotBCellMergedGR, file="tCellNotBCellGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(tCellNotBCellMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/tCellNotBCell.bed")

bCellNotTCellRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                                  presentOR = "Bcell", 
                                                  notPresentOR = c("CD4Tcell", "CD8Tcell"))]
bCellNotTCellGR = regionStringToGR(bCellNotTCellRegions)
save(bCellNotTCellGR, file="bCellNotTCellNotMergedGR.RData")
bCellNotTCellMergedGR = reduce(bCellNotTCellGR, min.gapwidth=mergeGapLength)
save(bCellNotTCellMergedGR, file="bCellNotTCellGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(bCellNotTCellMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/bCellNotTCell.bed")

# CD4 vs CD8 T cells
cd4NotCD8TCellRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                                     presentOR = "CD4Tcell", 
                                                     notPresentOR = "CD8Tcell")]
cd4NotCD8TCellGR = regionStringToGR(cd4NotCD8TCellRegions)
save(cd4NotCD8TCellGR, file="cd4NotCD8TCellNotMergedGR.RData")
cd4NotCD8TCellMergedGR = reduce(cd4NotCD8TCellGR, min.gapwidth=mergeGapLength)
save(cd4NotCD8TCellMergedGR, file="cd4NotCD8TCellGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(cd4NotCD8TCellMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/cd4NotCD8TCell.bed")

cd8NotCD4TCellRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                                      presentOR = "CD8Tcell", 
                                                      notPresentOR = "CD4Tcell")]
cd8NotCD4TCellGR = regionStringToGR(cd8NotCD4TCellRegions)
save(cd8NotCD4TCellGR, file="cd8NotCD4TCellNotMergedGR.RData")
cd8NotCD4TCellMergedGR = reduce(cd8NotCD4TCellGR, min.gapwidth=mergeGapLength)
save(cd8NotCD4TCellMergedGR, file="cd8NotCD4TCellGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(cd8NotCD4TCellMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/cd8NotCD4TCell.bed")

# NK cell vs other differentiated lymphoid cell
nkNotOtherDiffLymphoidRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                                      presentOR = "NKcell", 
                                                      notPresentOR = c("CD4Tcell", "CD8Tcell", "Bcell"))]
nkNotOtherDiffLymphoidGR = regionStringToGR(nkNotOtherDiffLymphoidRegions)
save(nkNotOtherDiffLymphoidGR, file="nkNotOtherDiffLymphoidNotMergedGR.RData")
nkNotOtherDiffLymphoidMergedGR = reduce(nkNotOtherDiffLymphoidGR, min.gapwidth=mergeGapLength)
save(nkNotOtherDiffLymphoidMergedGR, file="nkNotOtherDiffLymphoidGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(nkNotOtherDiffLymphoidMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/nkNotOtherDiffLymphoid.bed")

otherDiffLymphoidNotNKRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                                              presentOR = c("CD4Tcell", "CD8Tcell", "Bcell"), 
                                                              notPresentOR = "NKcell")]
otherDiffLymphoidNotNKGR = regionStringToGR(otherDiffLymphoidNotNKRegions)
save(otherDiffLymphoidNotNKGR, file="otherDiffLymphoidNotNKNotMergedGR.RData")
otherDiffLymphoidNotNKMergedGR = reduce(otherDiffLymphoidNotNKGR, min.gapwidth=mergeGapLength)
save(otherDiffLymphoidNotNKMergedGR, file="otherDiffLymphoidNotNKGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(otherDiffLymphoidNotNKMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/otherDiffLymphoidNotNK.bed")


# GMP vs MEP
gmpNotMEPRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                                      presentOR = "GMP", 
                                                      notPresentOR = "MEP")]
gmpNotMEPGR = regionStringToGR(gmpNotMEPRegions)
save(gmpNotMEPGR, file="gmpNotMEPNotMergedGR.RData")
gmpNotMEPMergedGR = reduce(gmpNotMEPGR, min.gapwidth=mergeGapLength)
save(gmpNotMEPMergedGR, file="gmpNotMEPGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(gmpNotMEPMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/gmpNotMEP.bed")

mepNotGMPRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                                 presentOR = "MEP", 
                                                 notPresentOR = "GMP")]
mepNotGMPGR = regionStringToGR(mepNotGMPRegions)
save(mepNotGMPGR, file="mepNotGMPNotMergedGR.RData")
mepNotGMPMergedGR = reduce(mepNotGMPGR, min.gapwidth=mergeGapLength)
save(mepNotGMPMergedGR, file="mepNotGMPGR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(mepNotGMPMergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/mepNotGMP.bed")

# present in common myeloid progenitor or common lymphoid progenitor but not in
# HSCs or terminally differented cells (picking CMP and CLP is somewhat arbitrary
# since there are other cell types in the middle of the hematopoiesis tree)
middleOfTree = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                             presentOR = c("CMP", "CLP"), 
                                             notPresentOR = c("HSC", "Bcell", "CD4Tcell", "CD8Tcell", "NKcell", "Mono", "Ery"))]
CMP_CLP_Not_HSC_Not_Diff_GR = regionStringToGR(middleOfTree)
save(CMP_CLP_Not_HSC_Not_Diff_GR, file="CMP_CLP_Not_HSC_Not_Diff_NotMergedGR.RData")
CMP_CLP_Not_HSC_Not_Diff_MergedGR = reduce(CMP_CLP_Not_HSC_Not_Diff_GR, min.gapwidth=mergeGapLength)
save(CMP_CLP_Not_HSC_Not_Diff_MergedGR, file="CMP_CLP_Not_HSC_Not_Diff_GR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(CMP_CLP_Not_HSC_Not_Diff_MergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/CMP_CLP_Not_HSC_Not_Diff.bed")

notMiddleOfTree = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat=openCellMat, 
                                             presentOR = c("HSC", "Bcell", "CD4Tcell", "CD8Tcell", "NKcell", "Mono", "Ery"), 
                                             notPresentOR = c("CMP", "CLP"))]
hsc_diff_Not_CMP_Not_CLP_GR = regionStringToGR(notMiddleOfTree)
save(hsc_diff_Not_CMP_Not_CLP_GR, file="hsc_diff_Not_CMP_Not_CLP_NotMergedGR.RData")
hsc_diff_Not_CMP_Not_CLP_MergedGR = reduce(hsc_diff_Not_CMP_Not_CLP_GR, min.gapwidth=mergeGapLength)
save(hsc_diff_Not_CMP_Not_CLP_MergedGR, file="hsc_diff_Not_CMP_Not_CLP_GR.RData")
RGenomeUtils::writeBed(df = MIRA:::grToDt(hsc_diff_Not_CMP_Not_CLP_MergedGR), filename = "hg19/hema_ATAC_GSE75384/regions/CMP_CLP_Not_HSC_Not_Diff.bed")



# Present in one specific cell type but in no others
cellType = colnames(openCellMat)
cellSpecificGRList = list()
noRegions = rep(TRUE, length(cellType))
regionSetName = paste0(cellType, "_specific")
for (i in seq_along(cellType)) {
    # present in one, not present in all others
    cellSpecRegions = regionString[subtypeRegionSet(notPresentMat=lowCellMat, presentMat = openCellMat, presentOR = cellType[i], notPresentOR = cellType[-i])]
    cellSpecificGRList[[i]] = regionStringToGR(cellSpecRegions)
    if (length(cellSpecificGRList[[i]]) > 0) {
        noRegions[[i]] = FALSE
        assign(paste0(regionSetName[i], "_NotMerged"), regionStringToGR(cellSpecRegions))
        save(list=paste0(regionSetName[i], "_NotMerged"), file=paste0(regionSetName[i],"_NotMerged.RData"))
        assign(regionSetName[i], reduce(get(paste0(regionSetName[i], "_NotMerged")), min.gapwidth=mergeGapLength))
        save(list=regionSetName[i], file=paste0(regionSetName[i], ".RData"))
        RGenomeUtils::writeBed(df = MIRA:::grToDt(get(paste0(regionSetName[i]))), 
                               filename = paste0("hg19/hema_ATAC_GSE75384/regions/", regionSetName[i], ".bed"))
    }
}
cellSpecificGRList = cellSpecificGRList[!noRegions]
names(cellSpecificGRList) = paste0(cellType[!noRegions], "_specific")
cellTypeSpecificGRList = GRangesList(cellSpecificGRList)
# sapply(X = cellTypeSpecificGRList, FUN = function(x) length(x))
save(cellTypeSpecificGRList, file= "cellTypeSpecificNotMergedGRList.RData")
cellTypeSpecificMergedGRList = lapply(X = cellTypeSpecificGRList, 
                                      FUN = function(x) reduce(x, min.gapwidth=mergeGapLength))
cellTypeSpecificMergedGRList = GRangesList(cellTypeSpecificMergedGRList)
save(cellTypeSpecificMergedGRList, file="cellTypeSpecificGRList.RData")

## liftover to hg38 genome with command line tool
# cd $RESOURCES/regions/LOLAExt/hg19/hematopoietic_ATACseq_GSE75384/regions/
# for FILE in $(ls)
# do
#   liftOver ./$FILE $RESOURCES/liftover/hg19ToHg38.over.chain $RESOURCES/regions/LOLAExt/hg38/hematopoietic_ATACseq_GSE75384/regions/$FILE $RESOURCES/regions/LOLAExt/hg38/hematopoietic_ATACseq_GSE75384/unmapped.txt
# done

# for the above region sets, merge adjacent regions

# allRegions = regionStringToGR(regionString)
# OL = findOverlaps(allRegions, allRegions)
# # checking whether there are any close regions which should perhaps be merged
# # originally each region is 500 bp
# allRegions2 = resize(x=allRegions, width = 1000, fix = "center")
# OL2 = findOverlaps(allRegions2, allRegions2)
# length(subjectHits(OL2)) == length(subjectHits(OL))

# # test out reduce function
# b = reduce(monoNotEryGR, min.gapwidth=501L) # gap of 200 is merged
# length(b)
# length(reduce(allRegions2, min.gapwidth=1L))
# test = GRanges(seqnames = rep("chr1", 4), ranges = IRanges(start = c(1, 21, 41, 61), end=c(20, 30, 59, 70)))
# reduce(test, min.gapwidth=1L)

# merge regions that have X bp or less separating them (only if both
# regions are open in a given cell type, otherwise don't merge: do this 
# as final step after getting custom region sets)

################################################################################
# visualize atac signal in defined cell type specific regions
# to confirm visually that they really are only open in one cell type

# normalized count distributions for each cell type
multiNiceHist(file = ffamlPlots(paste0(plotSubdir, "normalizedATACByCellType.pdf")), 
              dataDF = as.data.frame(mainCellMat), 
              colsToPlot = colnames(mainCellMat), xLabels = "Normalized ATAC counts", 
              binwidth = 1, boundary = 0, plotTitles = colnames(mainCellMat))

devtools::load_all(paste0(Sys.getenv("CODE"), "COCOA/"))
load("cellTypeSpecificGRList.RData")

# get average atac signal per region
byRegionCellType = lapply(cellTypeSpecificMergedGRList,FUN = function(x) averagePerRegion(signal = mainCellMat, 
                                                                       signalCoord = regionStringToGR(regionString), 
                                                                       regionSet = x, 
                                                                       signalCol = colnames(mainCellMat)))
byRegionCellType = lapply(X = byRegionCellType, as.data.frame)
byRegionCellType = lapply(X = byRegionCellType, function(x) x[, colnames(mainCellMat)])

nRSToPlot = length(cellTypeSpecificMergedGRList)
pdf(ffamlPlots(paste0(plotSubdir, "atacInCellTypeSpecificRegions.pdf")), width=11, height = 8.5 * nRSToPlot)
for (i in seq_along(byRegionCellType)) {
    multiHM = grid.grabExpr(draw(Heatmap(byRegionCellType[[i]], cluster_rows = TRUE, name = colnames(mainCellMat)[i])))
    pushViewport(viewport(y = unit((8.5 * nRSToPlot) - 
                                       (i - 1) * 8.5, "in"), 
                          height = unit(8, "in"), 
                          just = "top"))
    grid.draw(multiHM)
    popViewport()
}

dev.off()

# also look at the unnormalized signal in these regions
dim(regionCounts)
# because certain symbols in column names cause problems (and starting with # is problem)
colnames(regionCounts) <- paste0("sample_", gsub(pattern = "-", replacement = "_", x = colnames(regionCounts)))
byRegionSampleCounts = lapply(cellTypeSpecificMergedGRList,FUN = function(x) averagePerRegion(signal = regionCounts, 
                                                                                          signalCoord = regionStringToGR(regionString), 
                                                                                          regionSet = x, 
                                                                                          signalCol = colnames(regionCounts)))
byRegionSampleCounts = lapply(X = byRegionSampleCounts, as.data.frame)
byRegionSampleCounts = lapply(X = byRegionSampleCounts, function(x) x[, colnames(regionCounts)])

nRSToPlot = length(cellTypeSpecificMergedGRList)
pdf(ffamlPlots(paste0(plotSubdir, "atacInCellTypeSpecificRegionsBySample.pdf")), width=11, height = 20 * nRSToPlot)
for (i in seq_along(byRegionSampleCounts)) {
    multiHM = grid.grabExpr(draw(Heatmap(t(byRegionSampleCounts[[i]]), cluster_rows = TRUE, name = names(byRegionSampleCounts)[i])))
    pushViewport(viewport(y = unit((20 * nRSToPlot) - 
                                       (i - 1) * 20, "in"), 
                          height = unit(19.5, "in"), 
                          just = "top"))
    grid.draw(multiHM)
    popViewport()
}

dev.off()

#################################################################################
############################################################################
# a less supervised way of getting cell type specific regions



# first take out special conditions that are not counted as normal cell types
#' # Not in Figure1: pHSC (preleukemic HSC), Blast (leukemic blast), 
#' # LSC (leukemia stem cell), CD34 Cord Blood, (CD34 Bone Marrow?)
#' # for clustering analysis, taking out CD34 Cord Blood and Bone Marrow (low #
#' # and not a single cell type)
#' # also taking out leukemia related sets: pHSC, Blast, LSC
mainCellCounts = cellTypeCounts[, !(colnames(cellTypeCounts) %in% c("CD34 Cord Blood", 
                                                                    "CD34 Bone Marrow", "LSC", "Blast", "pHSC"))]

if (!all(cellType == colnames(mainCellCounts))) {
    warning()
}

cellType = colnames(mainCellCounts)
ncol(mainCellCounts)

pdf(file = paste0(Sys.getenv("PLOTS"), plotSubdir, "cellTypeCounts_hist_byCell.pdf"))
for (i in seq_along(cellType)) {
    hist(mainCellCounts[, cellType[i]], xlab = "Count", main=cellType[i]) 
}
dev.off()
apply(X = mainCellCounts, 2, max)
# a matrix defining ideal distribution, open in one cell type, closed in others
# columns of cellRefMat are the cell types, each row is a reference vector
cellRefMat = diag(length(cellType))
cellRefMat = scale(cellRefMat, scale=FALSE)
colnames(cellRefMat) <- cellType


# each row corresponds to a region
# for (i in 1:nrow(mainCellCounts)) {
#     # each reference vector corresponds to a column in cellTypeCor
#     for (j in 1:nrow(cellRefMat)) {
#         cellTypeCor[i, j] = cor(mainCellCounts[i, ], cellRefMat[j, ])
#     }
# }

# otherwise distributions with larger values and variance will have higher correlation
mainCellCounts = scale(mainCellCounts)

# same thing as above but with apply statements
simpleCache("cellTypeCor", {
    cellTypeCor = matrix(nrow=nrow(mainCellCounts), ncol=nrow(cellRefMat))
    colnames(cellTypeCor) = cellType
    for (i in 1:nrow(cellRefMat)) {
        cellTypeCor[, i] = apply(X = mainCellCounts, MARGIN = 1, FUN = function(x) cor(x, cellRefMat[i, ]) )
    }
    cellTypeCor
}, recreate=TRUE)
ComplexHeatmap::Heatmap(cellTypeCor[1:50000, ], cluster_columns = FALSE, cluster_rows=FALSE)
sum(colSums(cellTypeCor))


#################### Visualize
pdf(file = paste0(Sys.getenv("PLOTS"), plotSubdir, "cellTypeCor_hist.pdf"))
hist(cellTypeCor)
dev.off()

pdf(file = paste0(Sys.getenv("PLOTS"), plotSubdir, "cellTypeCor_hist_byCell.pdf"))
for (i in seq_along(cellType)) {
    hist(cellTypeCor[, cellType[i]], xlab = "Correlation", main=cellType[i], ylim = c(0, 200000)) 
}
dev.off()

# visualize the effect of different correlation cutoffs on how many cell type specific open 
# chromatin regions there are per cell
corCutoff = .7
cellCorBin = t(apply(X = cellTypeCor, MARGIN = 1, function(x) as.numeric(x > corCutoff)))
# cellCorBin = t(apply(X = cellTypeCor, MARGIN = 1, function(x) as.numeric(x < corCutoff & x > (corCutoff - 0.20))))
# hist(cellCorBin[ ,18])$counts
cellTypeF = factor(cellType, levels=cellType)
pdf(file = paste0(Sys.getenv("PLOTS"), plotSubdir, "cellTypeCorCutoff.pdf"))
cutoffVec = c(.5, .6, .7, .8, .9)
for (i in seq_along(cutoffVec)) {
    corCutoff = cutoffVec[i]
    cellCorBin = t(apply(X = cellTypeCor, MARGIN = 1, function(x) as.numeric(x > corCutoff)))
    numCellsOpen = rowSums(cellCorBin)
    hist(numCellsOpen, main = paste0("Cutoff=", cutoffVec[i]))
    plot(cellTypeF, colSums(cellCorBin), 
         main = paste0("Cutoff=", cutoffVec[i]), 
         xlab = "Cell Type", 
         ylab="Number of Cell Type Specific Regions")
    legend("topleft", cellType)
    
}
dev.off()
numCellsOpen = rowSums(cellCorBin)
hist(numCellsOpen)$counts
colSums(cellCorBin)

closeTogether = rep(0, length(cutoffVec))
closeTogetherP = rep(0, length(cutoffVec))
rangeVec = c(.1, .2, .3, .4)
closeTogetherDF = data.frame("0.1"=closeTogether, "0.2"=closeTogether, "0.3"=closeTogether, "0.4"=closeTogether)
closeTogetherPDF = data.frame("0.1"=closeTogether, "0.2"=closeTogether, "0.3"=closeTogether, "0.4"=closeTogether)
# pdf(file = paste0(Sys.getenv("PLOTS"), plotSubdir, "cellTypeCor_Cutoff_robustness.pdf"))
for (iOuter in seq_along(rangeVec)) {
    
    for (i in seq_along(cutoffVec)) {
        corCutoff = cutoffVec[i]
        cellCorBin = t(apply(X = cellTypeCor, MARGIN = 1, function(x) as.numeric(x < corCutoff & x > (corCutoff - rangeVec[iOuter]))))
        numCellsOpen = rowSums(cellCorBin)
        closeTogether[i] = sum(numCellsOpen > 1)
        # proportion
        cellCorBin = t(apply(X = cellTypeCor, MARGIN = 1, function(x) as.numeric(x > corCutoff & x > (corCutoff - rangeVec[iOuter]))))
        numCellsOpen = rowSums(cellCorBin)
        closeTogetherP[i] = sum(numCellsOpen > 1) / sum(numCellsOpen > 0)
    }
    closeTogetherDF[, iOuter] = closeTogether
    closeTogetherPDF[, iOuter] = closeTogetherP
}
closeTogetherDF = cbind(closeTogetherDF, corThreshold=cutoffVec)
closeTogetherDF = tidyr::gather(data = closeTogetherDF, key=rangeThreshold, value=numRegions, -corThreshold)

threshRobustPlot = ggplot(data = closeTogetherDF, aes(x = corThreshold, y=numRegions)) + geom_line(aes(col=rangeThreshold))
ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "cellTypeCor_Cutoff_robustness.pdf"), plot = threshRobustPlot)

# are most of the regions normally distributed in terms of correlation values?
pdf(file=paste0(Sys.getenv("PLOTS"), plotSubdir, "cellTypeCorCutoff.pdf"))
for (i in 1:1000) {
    hist(cellTypeCor[i, ], main = "Hist of individual region",
         xlab = "Correlation", 
         ylab="Number of Cell Types")
}
dev.off()

###########################
# using two threshold: must be above a certain correlation and the next highest cell type
# must below the correlation of the highest by a certain amount
corThresh = 0.5
diffThresh = 0.3 # difference of at least x
maxRegionCor = apply(X = cellTypeCor, 1, max) # cellTypeCor[, "max"]
secondHighest = apply(X = cellTypeCor, 1, function(x) sort(x, decreasing=TRUE)[2])
aboveDiffVec = (maxRegionCor - secondHighest) >= diffThresh
sum(aboveDiffVec)

cellCorBin = cellTypeCor > corThresh
# hist(rowSums(cellCorBin))
# what we define as cell type specific regions
cellCorBin = cellCorBin & aboveDiffVec
colSums(cellCorBin)
# plot(colSums(cellCorBin))
# get regions for each cell type (column)
hist(cellTypeCor[cellCorBin[, 13], 13])

# plot(cutoffVec, closeTogether, 
#      ylab = "Number of regions with at least 2 correlations within 0.2 of each other", main="Cutoff Robustness", xlab="Threshold")
# dev.off()
cellTypeSpecificGRList = GRangesList()
# make region sets
for (i in 1:ncol(cellCorBin)) {
    cellTypeSpecificGRList[[i]] = regionStringToGR(regionString[cellCorBin[, i]])
}
names(cellTypeSpecificGRList) = colnames(cellCorBin)

testInd = grep(pattern = "chr1_2076559_2077059", x = regionString)
cellCorBin[testInd, ]
cellTypeCor[testInd, ]

# also calculate correlation between cell type vector for more than one cell type and region
# columns of cellRefMat are the cell types, each row is a reference vector
# for all combinations up to ??
# all combinations of 2
cellRefMat = cellRefMat

# all combinations of 3

# all combinations of 4

##############
corThresh = 0.5
diffThresh = .3
# cell type specific closed chromatin
minRegionCor = apply(X = cellTypeCor, 1, min) # cellTypeCor[, "max"]
secondLowest = apply(X = cellTypeCor, 1, function(x) sort(x, decreasing=FALSE)[2])
belowDiffVec = (secondLowest - minRegionCor) >= diffThresh
sum(belowDiffVec)

cellCorBin = cellTypeCor < -corThresh
# hist(rowSums(cellCorBin))
# what we define as cell type specific regions
cellCorBin = cellCorBin & belowDiffVec
colSums(cellCorBin)



#################### Notes on the dataset
#' # Not in Figure1: pHSC (preleukemic HSC), Blast (leukemic blast), 
#' # LSC (leukemia stem cell), CD34 Cord Blood, (CD34 Bone Marrow?)
#' # for clustering analysis, taking out CD34 Cord Blood and Bone Marrow (low #
#' # and not a single cell type)
#' # also taking out leukemia related sets: pHSC, Blast, LSC
#' mainCellMat = cellTypeCounts[, !colnames(cellTypeCounts) %in% c("CD34 Cord Blood", 
#'                                   "CD34 Bone Marrow", "LSC", "Blast", "pHSC")]
#' 
#' # normalizing since each cell type has a different mean/median count (column normalization)
#' mainCellMat = scale(mainCellMat, center = FALSE, scale = apply(mainCellMat, 2, max))
#' 
#' # doing k means clustering (columns are features (cell type))
#' # no hierarchical clustering (distance matrix is too big)
#' set.seed(10)
#' regionKClust1 = kmeans(mainCellMat, centers = 100, iter.max = 1000, algorithm = "MacQueen")
#' # # checking how well clustering worked
#' # # regionDist = dist(mainCellMat) # too much memory required 
#' # # cluster::silhouette(regionKClust$cluster, regionDist)
#' # set.seed(11)
#' # regionKClust2 = kmeans(mainCellMat, centers = 75, iter.max = 1000, algorithm = "MacQueen")
#' # set.seed(100)
#' # regionKClust3 = kmeans(mainCellMat, centers = 50, iter.max = 1000, algorithm="MacQueen")
#' # regionKClust4 = kmeans(mainCellMat, centers = 35, iter.max = 1000, algorithm="MacQueen")
#' 
#' # # visualize clusters
#' # a = as.data.frame(regionKClust1$centers)
#' # a = cbind(a, cluster_num=1:13)
#' # # convert to long format first 
#' # b = gather(a, cell_type, value, HSC:Ery, factor_key = TRUE)
#' 
#' #' This function takes the object of cluster averages (centroid) and
#' #' makes a separate bar plot for each cluster with the dimensions on the x axis
#' #' and the value on the y axis
#' 
#' #' @param df data.frame with column for cluster number, dimensions in which cluster was done and 
#' #' value in that dimension (long format) eg 1, dim1, 5
#' #' @param xCol character object, name of column for x axis
#' #' @param yCol character object, name of column for y axis
#' #' @param clusterCol name of column with cluster number
#' #' @param cluster_num only if facetWrap =FALSE, a single cluster to plot
#' #' to input matrix directly from the kmeans return object set clusterCol to NULL
#' #' @examples plotClusterDims(regionKClust1$centers, "cell_type", "value", clusterCol = NULL, facetWrap=TRUE)
#' #' plotClusterDims(regionKClust1$centers, "cell_type", "value", clusterCol = NULL, facetWrap=FALSE, cluster_num = 20)
#' #' @export
#' plotClusterDims = function(df, xCol, yCol, clusterCol, cluster_num = NULL, facetWrap=TRUE) {
#'     if (is.null(clusterCol)) {
#'         df = as.data.frame(df)
#'         df = cbind(df, clusterNum = 1:nrow(df))
#'         df = gather(df, tmp1, tmp2, colnames(df)[1:(ncol(df) - 1)], factor_key = TRUE)
#'         clusterCol = "clusterNum"
#'         names(df) = c("clusterNum", xCol, yCol)
#'     }
#'     if (facetWrap) {
#'         ggplot(data=df, aes(x=get(xCol), y=get(yCol))) + geom_bar(stat = "identity") + facet_wrap(clusterCol)
#'     } else {
#'         ggplot(data=dplyr::filter(df, get(clusterCol) == cluster_num), aes(x=get(xCol), y=get(yCol))) + geom_bar(stat = "identity") #+ facet_wrap("cluster_num") 
#'     }
#' }
#' 
#' pK1 = plotClusterDims(regionKClust1$centers, "cell_type", "value", clusterCol = NULL, facetWrap=TRUE)
#' pK1
#' plotClusterDims(regionKClust2$centers, "cell_type", "value", clusterCol = NULL, facetWrap=TRUE)
#' plotClusterDims(regionKClust3$centers, "cell_type", "value", clusterCol = NULL, facetWrap=TRUE)
#' plotClusterDims(regionKClust4$centers, "cell_type", "value", clusterCol = NULL, facetWrap=TRUE)
#' plotClusterDims(regionKClust1$centers, "cell_type", "value", clusterCol = NULL, facetWrap=FALSE, cluster_num = 65)
#' table(regionKClust1$cluster)
#' ggsave(pK1, filename = paste0(Sys.getenv("PLOTS"), "hemaClust100.pdf"), device = "pdf")
#' 
#' # determine which k is best, factoextra package
#' #fviz_nbclust(x = mainCellMat, FUNcluster = kmeans, method = "wss")+
#'     # geom_vline(xintercept = 4, linetype = 2)+
#'     # labs(subtitle = "Elbow method")
#' # fails because it tries to make a dist matrix, not enough memory
#' 
#' # manually choose clusters to get region sets for
#' # choosing clusters: 10 (CLP, 42 regions), 86 (MEP, CMP, 116 regions), 60 (1940 regions), 2 (2150 regions), 43 (444 regions), 5 (924 regions)
#' # get regions for those clusters and convert each to a GRanges region set
#' clustRegionList = list()
#' chosenClust = c(2, 5, 10, 43, 60, 86)
#' for (i in seq_along(chosenClust)) {
#'     # regions for this cluster
#'     theseRegions = regionString[regionKClust1$cluster == chosenClust[i]]
#'     chr = sub(pattern = "_.+_.+", replacement = "", theseRegions)
#'     start = sub(pattern = "chr.{1,2}_", "", theseRegions)
#'     start = as.numeric(sub(pattern = "_.+$", "", start))
#'     end = as.numeric(sub(pattern = ".+_.+_", "", theseRegions))
#'     clustRegionList[[i]] = data.frame(chr, start, end)
#' }
#' names(clustRegionList) = as.character(chosenClust)
#' 
#' clustRegionGRList = lapply(X = clustRegionList, MIRA:::dtToGr)
#' 
#' # GRangesList with a separate region set for each sample
#' # top regions but only in respect to the cell type itself, not to other cell types
#' save(clustRegionGRList, file="GSE74912_Test_Clusters_Region_Sets.RData")
