# permutation test by shuffling sample labels

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "19-permTestBRCA_ATAC_Protein"
plotSubdir = "19-permBRCA_ATAC_Protein/"

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 250

######################################################################
# required inputs to permutation test


# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg38")

colsToAnnotate = #protein columns

dataID = "brcaATACProtein"

############################################################################

source(ffProjCode("src/runPermTest.R"))