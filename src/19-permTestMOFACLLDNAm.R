# permutation test by shuffling sample labels

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
devtools::load_all(ffCode("COCOA/"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "19-permMOFACLLDNAm"
plotSubdir = "19-permMOFACLLDNAm/"
dataID = "CLL196MOFA"

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 300

######################################################################
variationMetric = "cov"
dataID = paste0(dataID, "_", variationMetric)

# assigns methylMat, signalCoord, latentFactorMat to global environment
loadMOFAData()
genomicSignal = methylMat
sampleLabels = data.frame(latentFactorMat)

# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg19")
simpleCache(paste0("rsScore_", dataID), assignToVariable = "realRSScores")
row.names(realRSScores) = realRSScores$rsName
sharedRSNames = names(GRList)[names(GRList) %in% realRSScores$rsName]
GRList = GRList[sharedRSNames]
realRSScores = realRSScores[sharedRSNames, ]


# colsToAnnotate = paste0("LF", c(1:3, 5:7, 9))
colsToAnnotate = paste0("LF", 1:10)



#####################################################################
# permutation test for significance
# requires: nPerm, sampleLabels, genomicSignal, signalCoord, GRList, colsToAnnotate
# dataID, variationMetric (optional, default="cor")
source(ffProjCode("src/runPermTest.R"))

####################################################################

# multiNiceHist(file = ffPlot(paste0(plotSubdir, "nullDistMOFACorPermRS2250.pdf")), dataDF = nullDistList[[2250]], 
#               colsToPlot = colsToAnnotate, xLabels = "COCOA score (absolute correlation)", 
#               binwidth = 0.001, yLabel = "Number of permutation results", 
#               plotTitles = paste0("Null distribution of COCOA scores, ", colsToAnnotate),
#               ggExpr = "+xlim(0,0.15)")
# 
# hist(nullDistList[[2260]]$LF1)
# 
# simpleCache(paste0("mofaPermPVals", dataID), {
#     rsPVals = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
#                           calcCols=colsToAnnotate, whichMetric = "pval")
#     rsPVals
# })
# multiNiceHist(file = ffPlot(paste0(plotSubdir, "pValDistMOFACorPerm.pdf")), dataDF = rsPVals, 
#               colsToPlot = colsToAnnotate, xLabels = "p-value", 
#               binwidth = 0.005, yLabel = "Number of region sets", 
#               plotTitles = paste0("Distribution of region set p-values, ", colsToAnnotate),
#               ggExpr = "+xlim(0,1)+ylim(0, 2270)")
# 
# 
# View(rsPVals[which(rsPVals$LF5 < 0.01), ])
# 
# simpleCache(paste0("mofaPermZScores", dataID), { 
#     rsZScores = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
#                             calcCols=colsToAnnotate, whichMetric = "zscore")
#     rsZScores
#     
# })
# View(rsZScores[which(rsZScores$LF2 > 7), ])
# 
# multiNiceHist(file = ffPlot(paste0(plotSubdir, "zScoreDistMOFACorPerm.pdf")), dataDF = rsZScores, 
#               colsToPlot = colsToAnnotate, xLabels = "z score", 
#               binwidth = 1, yLabel = "Number of region sets", 
#               plotTitles = paste0("Distribution of region set z scores, ", colsToAnnotate),
#               ggExpr = "+xlim(-3,40)")
# 
# hist(rsZScores$LF5)
# 
# topRSInd = rsRankingIndex(rsScores = rsZScores, PCsToAnnotate = colsToAnnotate)
# 
# topRSAnnoList = list()
# for (i in seq_along(colsToAnnotate)) {
#     topRSAnnoList[[i]] = data.frame(rsName=rsZScores$rsName[topRSInd[1:20, colsToAnnotate[i]]], 
#                                     rsDescription=rsZScores$rsDescription[topRSInd[1:20, colsToAnnotate[i]]])
#     names(topRSAnnoList[[i]]) <- paste0(names(topRSAnnoList[[i]]), "_", colsToAnnotate[i])
# }
# 
# write.csv(topRSAnnoList, file = ffProc(paste0("COCOA_paper/analysis/sheets/topRSMOFACorPerm.csv")))
# for (i in 2:length(topRSAnnoList)) {
#     write.csv(topRSAnnoList[[i]], file = ffProc(paste0("COCOA_paper/analysis/sheets/topRSMOFACorPerm.csv")), append = TRUE)
#     
# }
