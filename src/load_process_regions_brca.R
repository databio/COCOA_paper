
library(LOLA)
library(GenomicRanges)
library(data.table)
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-genericFunctions.R"))

# reading in the region sets
# load LOLA database
lolaPath1 = paste0(Sys.getenv("RESOURCES"), "/regions/LOLACore/hg38/")
regionSetDB = loadRegionDB(lolaPath1)
loRegionAnno = regionSetDB$regionAnno
lolaCoreRegionAnno = loRegionAnno
# a549Ind = grep("a549", loRegionAnno$cellType, ignore.case = TRUE)
sheff_dnaseInd = grep("sheffield_dnase", loRegionAnno$collection, ignore.case = TRUE)
# mcf7Ind = grep("mcf-7", loRegionAnno$cellType, ignore.case = TRUE)
# k562Ind = grep("k562", loRegionAnno$cellType,  ignore.case = TRUE)
# GRList = GRangesList(regionSetDB$regionGRL[c(a549Ind, mcf7Ind)])
fInd1 = filterFetal(lolaCoreRegionAnno)
lolaCoreRegionAnno = lolaCoreRegionAnno[-sort(unique(c(sheff_dnaseInd, fInd1)))]
GRList1 = GRangesList(regionSetDB$regionGRL[-sort(unique(c(sheff_dnaseInd, fInd1)))])
allRegionAnno = lolaCoreRegionAnno

# ROADMAP Epigenome project and Jaspar motifs
lolaPath2 = paste0(Sys.getenv("RESOURCES"), "/regions/LOLAExt/hg38/")
regionSetDB2 = loadRegionDB(lolaPath2, useCache = TRUE)
loRegionAnno2 = regionSetDB2$regionAnno
roadmapRegionAnno = loRegionAnno2[loRegionAnno2$collection == "roadmap_epigenomics", ]
GRList2 = GRangesList(regionSetDB2$regionGRL[loRegionAnno2$collection == "roadmap_epigenomics"])
fInd2 = filterFetal(roadmapRegionAnno)
roadmapRegionAnno = roadmapRegionAnno[-fInd2]
GRList2 = GRList2[-fInd2]
allRegionAnno = rbind(allRegionAnno, roadmapRegionAnno)

# processing Jaspar motif regions
# resizing regions and filtering to only open chromatin in K562
# size was ~999
motifRegionAnno = loRegionAnno2[loRegionAnno2$collection == "jaspar_motifs", ]
allRegionAnno = rbind(allRegionAnno, motifRegionAnno)
GRList3 = GRangesList(regionSetDB2$regionGRL[loRegionAnno2$collection == "jaspar_motifs"])
# filtering based on chromatin accessibility data from blood/AML
# load(paste0(Sys.getenv("PROCESSED"), ))
# GRList3 = getOLRegions(GRList = GRList3, intGR=bloodAccessibleRegions, removeOL = FALSE)
# GRList3 = GRangesList(GRList3)
#I'm not sure that center is where motif is: GRList3 = resize(GRList3, width = 200, fix="center") 

hemaATACRegionAnno = loRegionAnno2[loRegionAnno2$collection == "hematopoietic_ATACseq_GSE75384", ]
GRList4 = GRangesList(regionSetDB2$regionGRL[loRegionAnno2$collection == "hematopoietic_ATACseq_GSE75384"])
allRegionAnno = rbind(allRegionAnno, hemaATACRegionAnno)

# combine into one GRList
GRList = c(GRList1, GRList2, GRList3, GRList4)

#adding ER Chipseq dataset
# erSet = fread(paste0(Sys.getenv("CODE"), "COCOA_paper/inst/extdata/",
#                      "GSM2305313_MCF7_E2_peaks_hg38.bed"))
# setnames(erSet, c("V1", "V2", "V3"), c("chr", "start", "end"))
# GRList = c(GRangesList(MIRA:::dtToGr(erSet)), GRList)
# names(GRList) = c("GSM2305313_MCF7_E2_peaks_hg38.bed", names(GRList))
#######################################################################
# c("GSM2305313_MCF7_E2_peaks_hg38.bed"
rsName = c(lolaCoreRegionAnno$filename,
           roadmapRegionAnno$filename,
           motifRegionAnno$filename,
           hemaATACRegionAnno$filename)
# c("ER ChIPseq, MCF7 cell line, estradiol stimulation"
rsDescription = c(lolaCoreRegionAnno$description,
                  roadmapRegionAnno$description,
                  motifRegionAnno$description,
                  hemaATACRegionAnno$description)
rsCollection = c(lolaCoreRegionAnno$collection, roadmapRegionAnno$collection, 
                 motifRegionAnno$collection, hemaATACRegionAnno$collection)

#################################################################
# cleaning up since there were many large objects
if (exists("GRList")) {
    rm(list = c("GRList1", "GRList2", "GRList3", "regionSetDB", "regionSetDB2"))
    gc()
} else {
    warning("Loading GRList was not successful.")
}
