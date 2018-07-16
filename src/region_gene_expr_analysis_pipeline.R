# Gene expression related to regions that vary in DNA methylation 



# Link regions of interest with genes
# various options for this
# 1. Elmer::GetNearGenes (closest x genes)
# 2. genomation (overlap with known gene annotation, different parts of gene)
# ?genomation::getAssociationWithTSS, ?genomation::readTranscriptFeatures 
# 3. Genehancer (known enhancer gene links)
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=hub_148281&hgta_track=hub_148281_geneHancer_hg19&hgta_table=hub_148281_geneHancer_hg19&hgta_doSchema=describe+table+schema
# https://genecards.weizmann.ac.il/geneloc/gh_hub/hg19/ucsc_genehancer_hg19.bb
# UCSC program bigBedToBed 

# Does expression of genes vary with methylation and also with corresponding PC?


# Gene set enrichment analysis for genes that vary along with these regions?

# connecting region sets to nearest genes to understand significance of region sets

############################################################################33
# old code 


# get top 10 region sets from each PC
topRSInd = unique(unlist(rsEnSortedInd[1:3, ]))
regionSetList = GRList[topRSInd]

# only regions that covered at least one measured cytosine
# top region sets
coveredRSList = getOLRegions(GRList = regionSetList, intGR= MIRA:::dtToGr(bigSharedC$coordinates))
# all region sets
covAllRSList = getOLRegions(GRList = GRList, intGR= MIRA:::dtToGr(bigSharedC$coordinates))

library(GenomicFeatures) # bioc help page 76512
hg19.refseq.db <- makeTxDbFromUCSC(genome="hg19", tablename="knownGene") # could use "knownGene"
# supportedUCSCtables
refseq.genes <- genes(hg19.refseq.db)
# refseq.genes<- transcripts(hg19.refseq.db)
t3 = GenomicRanges::nearest(x = coveredRSList[[10]], subject = refseq.genes)


# might need to use apply
test = getNearest(what = "gene", probes= coveredRSList[[10]])
names(coveredRSList[[1]]) <- NULL
t2 = getNearestGene(probes = coveredRSList[[1]])

# only want genes related to highly variable regions, not all regions
# identify highly variable regions for each PC

# also screen based on closeness (linear genomic distance)?

######### gene expression and PCA score/DNAm 
# testing whether gene expression at related genes is correlated with 
# DNA methylation at regions in the region set or with PC score (for original
# PC or for PC from only CpG in region set)