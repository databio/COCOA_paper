# curate region set databases


loadGRList(genomeV = "hg19")

#############################################################################
# get an idea of the what transcription factors are in the database


# list of TFs downloaded from: http://fantom.gsc.riken.jp/5/sstar/Browse_Transcription_Factors_hg19
fantomFile = read.csv(paste0(Sys.getenv("RESOURCES"), "Browse Transcription Factors hg19 - resource_browser.csv"), 
                      stringsAsFactors = FALSE) 
tfNames = fantomFile$Symbol
tfAlias = rep("", length(tfNames))

# alternate names
tfAlias[grep(pattern = "SPI1", x = tfNames)] = "|PU1|PU.1"
tfAlias[grep(pattern = "TAL1", x = tfNames)] = "|SCL"

# includes all TFs in one search
allTFPattern = paste0(tfNames, collapse = "|")

length(grep(pattern = allTFPattern,x =  rsName, ignore.case = TRUE))
rsCollection == "roadmap_epigenomics"

# https://stackoverflow.com/questions/2969315/rhow-to-get-grep-to-return-the-match-rather-than-the-whole-string

tfNames = paste0(tfNames, tfAlias)
matches = list()
for (i in seq_along(tfNames)) {
    pat <- paste0(".*(", tfNames[i], ").*")
    matches[[i]] = sub(pattern = pat, replacement = "\\1", x = rsName[grepl(pat, rsName)], ignore.case = TRUE)
}
matches = unlist(matches)
sort(table(matches))
# parentheses in the regular expression seem to interfere with what 
# it natches (e.g. not matching "Ctcf")