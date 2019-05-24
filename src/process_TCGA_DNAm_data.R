library(data.table)
library(simpleCache)

######### process DNA methylation data ###################################
# normalizing names and folder structures after download with gdc client
setwd(paste0(Sys.getenv("DATA"),"/tcga/DNA_methylation/450k/"))

# gdc manifest has a key matching folder names with file names
gdc_man = fread("gdc_manifest.2018-03-14.txt")
tcgaID = regmatches(gdc_man$filename, regexpr(pattern = "TCGA-..-....", gdc_man$filename))
# get cancer project/type
tcgaProj = regmatches(gdc_man$filename, regexpr(pattern = "_.{2,6}\\.", gdc_man$filename))
tcgaProj = gsub(pattern = "_", replacement = "", tcgaProj)
tcgaProj = gsub(pattern = "\\.", replacement = "", tcgaProj)
# table(tcgaProj)
# some patients have multiple distinct metylation arrays
# sometimes they have a normal tissue array
if (length(tcgaID) != length(unique(tcgaID))) {
    warning("Some patients are represented more than once.")
    
    # figuring out which patients are represented more than once
    idt = table(tcgaID)
    # head(sort(idt, decreasing = TRUE))
    multiInd = which(idt > 1)
    multiPatient = names(multiInd)
    # which cancers have patients with multiple methylation files
    # table(tcgaProj[multiInd])
    # 31 out of 33 cancer types have some patients with multiple files
}

# creating annotation sheet that has file path, cancer type, patient ID
file_path = paste0(gdc_man$id, "/", gdc_man$filename)
cancer_type = tcgaProj
subject_ID = tcgaID
dnamAnno = data.table(subject_ID, cancer_type, file_path)
# restricting to just BRCA data for my subproject
# also excluding patients with multiple files
setCacheDir(cacheDir = paste0(Sys.getenv("PROCESSED"), "/brca_PCA/RCache/"))
uBRCAAnno = dnamAnno[ !(dnamAnno$subject_ID %in% multiPatient) & dnamAnno$cancer_type == "BRCA", ]
uBRCAAnno = uBRCAAnno[order(subject_ID), ]
brcaID = uBRCAAnno$subject_ID
setwd(paste0(Sys.getenv("DATA"), "/tcga/DNA_methylation/450k/https:/api.gdc.cancer.gov/data"))
# add pattern = \\.txt$
methylFiles = uBRCAAnno$file_path
#MIRA:::setLapplyAlias(cores= 2)
methylList = lapply(X = methylFiles, FUN = fread)
names(methylList) = brcaID


####################################################
#' format 450k methylation files
#' does operations by reference so methylList object will be altered
#' @param methylList 1 list item for each sample/data matrix.
#' @param colsToKeep character vector for columns to keep. 
#' allows option to keep extra metadata columns.
#' @param newColNames same length as colsToKeep. new names for those columns
#' in respective order.
#' @param roundMethyl boolean. round beta values to 2 numbers after the decimal 0.xx
formatTCGA450k = function(methylList, 
                                  colsToKeep=c("Beta_value", 
                                                "Chromosome",
                                                "Start"), 
                                  newColNames=c("methylProp", "chr", 
                                                "start"),
                          roundMethyl=TRUE){
    methylList = data.table::copy(methylList)
    defaultCols = c("Composite Element REF", 
                    "Beta_value",
                    "Chromosome",
                    "Start",
                    "End",
                    "Gene_Symbol", 
                    "Gene_Type", 
                    "Transcript_ID", 
                    "Position_to_TSS",
                    "CGI_Coordinate",
                    "Feature_Type")
    # if (colToKeep == "all") {
    #     colsToKeep = defaultCols
    # }
    colsToDiscard = setdiff(defaultCols, colsToKeep)
    # delete some columns 
    methylList = lapply(X = methylList, FUN = function(x) x[, (colsToDiscard) :=NULL])
    methylList = lapply(X = methylList, FUN = setnames, 
                                        old = colsToKeep, 
                                        new=newColNames)
    # methylList = lapply(X = methylList, FUN = setcolorder, c("chr", "start", "end", "methylProp"))
    # methylList = lapply(X = methylList, FUN = function(x) x[, coverage := 1])
    if (roundMethyl) {
        methylList = lapply(X = methylList, FUN = function(x) x[, methylProp := round(methylProp, 2)])
    }
    
    return(methylList[])
}

################## 

methylList = formatTCGA450k(methylList)
methylList = lapply(methylList, function(x) setorder(x, chr, start))
simpleCache("methylDTList", {
    methylList
})
coordinateDTAll = copy(methylList[[1]])[, list(chr, start)]
coordinateDT = copy(coordinateDTAll)
methylList = lapply(methylList, function(x) x[, c("chr", "start") := NULL])
methylMat = do.call(cbind, methylList)
setnames(methylMat, names(methylList)) 

# screening out XY chromosome and "*" 
chrScreenInd = coordinateDT$chr %in% c("*", "chrX", "chrY")
methylMat = methylMat[!chrScreenInd, ]
coordinateDT = coordinateDT[!chrScreenInd, ]
# screen out any rows with NA
naInd = apply(methylMat, 1, function(x) any(is.na(x)))
methylMat = methylMat[!naInd, ]
coordinateDT = coordinateDT[!naInd, ]
methylMat = as.matrix(methylMat)

simpleCache("combinedBRCAMethyl_noXY", {
    combinedBRCAMethyl_noXY = list(coordinates = coordinateDT,
                                   methylProp = methylMat)
    combinedBRCAMethyl_noXY
})


###########################################################################


