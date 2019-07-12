# permutation test by shuffling sample labels

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "19-applyPerm"
plotSubdir = "19-applyPerm/"

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 250

######################################################################

# assigns signalMat, signalCoord, loadingMat, pcScores
loadBRCADNAm()
genomicSignal = signalMat
sampleLabels = pcScores

# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg38")

colsToAnnotate = paste0("PC", 1:10)

dataID = "brcaDNAm657"

############################################################################



indList = list()
# generate random indices for shuffling of samples
for (i in 1:nPerm) {
    indList[[i]] = sample(1:nrow(sampleLabels), nrow(sampleLabels))
}

# # plug random indices into parallelized function
# rsPermScores = COCOA:::lapplyAlias(X= indList, FUN=function(x) corPerm(randomInd=x, 
#                                                         genomicSignal=methylMat, 
#                                                         signalCoord=signalCoord, 
#                                                         GRList=GRList, 
#                                                         calcCols=colsToAnnotate,
#                                                         sampleLabels=latentFactors))

rsPermScores = list()
for (i in seq_along(indList)) {
    
    rsPermScores[[i]] = corPerm(randomInd=indList[[i]], 
                                genomicSignal=genomicSignal, 
                                signalCoord=signalCoord, 
                                GRList=GRList, 
                                calcCols=colsToAnnotate,
                                sampleLabels=sampleLabels)
    message(i)
    if ((i %% 50) == 0) {
        save(rsPermScores, file = ffProc(paste0("COCOA_paper/RCache/rsPermScores_", dataID, ".RData")))
    }
    
}

save(rsPermScores, file = ffProc(paste0("COCOA_paper/RCache/rsPermScores_", dataID, ".RData")))

# # get score null distribution for each region set
# # one null distribution for each region set
# nullDistList = lapply(X = seq_along(GRList),
#                       FUN = function(x) extractNullDist(resultsList=rsPermScores, rsInd = x))
# 
# multiNiceHist(file = ffPlot(paste0(plotSubdir, "nullDistMOFACorPermRS2250.pdf")), dataDF = nullDistList[[2250]], 
#               colsToPlot = colsToAnnotate, xLabels = "COCOA score (absolute correlation)", 
#               binwidth = 0.001, yLabel = "Number of permutation results", 
#               plotTitles = paste0("Null distribution of COCOA scores, ", colsToAnnotate),
#               ggExpr = "+xlim(0,0.15)")
# 
# hist(nullDistList[[2260]]$LF1)
# 
# rsPVals = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
#                       calcCols=colsToAnnotate, whichMetric = "pval")
# multiNiceHist(file = ffPlot(paste0(plotSubdir, "pValDistMOFACorPerm.pdf")), dataDF = rsPVals, 
#               colsToPlot = colsToAnnotate, xLabels = "p-value", 
#               binwidth = 0.005, yLabel = "Number of region sets", 
#               plotTitles = paste0("Distribution of region set p-values, ", colsToAnnotate),
#               ggExpr = "+xlim(0,1)+ylim(0, 2270)")
# 
# 
# View(rsPVals[which(rsPVals$LF5 < 0.01), ])
# 
# rsZScores = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
#                         calcCols=colsToAnnotate, whichMetric = "zscore")
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
