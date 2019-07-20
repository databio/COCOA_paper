# permutation test by shuffling sample labels
# expected inputs: 
#' @param nPerm numeric. The number of permutations to do.
#' @param sampleLabels data.frame/matrix. Sample labels/values that 
#' you are running COCOA to find region sets associated with. These 
#' values will be shuffled for the permutation test. Rows are samples.
#' Each column is a sample label.
#' @param genomicSignal
#' @param signalCoord
#' @param GRList
#' @param colsToAnnotate character. The column names of `sampleLabels` that
#' you want to test.
#' @param dataID character. A unique identifier for this dataset 
#' (for saving results)


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
        save(rsPermScores, file = ffProc(paste0("COCOA_paper/RCache/rsPermScores_", 
                                                dataID, ".RData")))
    }
    
}

save(rsPermScores, file = ffProc(paste0("COCOA_paper/RCache/rsPermScores_", 
                                        dataID, ".RData")))



nullDistList = lapply(X = seq_along(GRList),
                      FUN = function(x) extractNullDist(resultsList=rsPermScores, rsInd = x))

# just an example of the null distributions for a single region set (arbitrarily rs1)
multiNiceHist(file = ffPlot(paste0(plotSubdir, "nullDistRS1", dataID, ".pdf")), dataDF = nullDistList[[1]], 
              colsToPlot = colsToAnnotate, xLabels = "COCOA score (absolute correlation)", 
              binwidth = 0.001, yLabel = "Number of permutation results", 
              plotTitles = paste0("Null distribution of COCOA scores, ", colsToAnnotate),
              ggExpr = "+xlim(0,0.25)")


simpleCache(paste0("permPVals", dataID), {
    rsPVals = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
                          calcCols=colsToAnnotate, whichMetric = "pval")
    rsPVals
}, recreate = TRUE, reload = TRUE)
multiNiceHist(file = ffPlot(paste0(plotSubdir, "pValDist", dataID, ".pdf")), dataDF = rsPVals, 
              colsToPlot = colsToAnnotate, xLabels = "p-value", 
              binwidth = 0.005, yLabel = "Number of region sets", 
              plotTitles = paste0("Distribution of region set p-values, ", colsToAnnotate),
              ggExpr = paste0("+xlim(0,1)+ylim(0, ", nrow(rsPVals), ")"))

simpleCache(paste0("permZScores", dataID), { 
    rsZScores = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
                            calcCols=colsToAnnotate, whichMetric = "zscore")
    rsZScores
    
}, recreate = TRUE, reload=TRUE)
# View(rsZScores[which(rsZScores$ > 7), ])

multiNiceHist(file = ffPlot(paste0(plotSubdir, "zScoreDist", dataID, ".pdf")), dataDF = rsZScores, 
              colsToPlot = colsToAnnotate, xLabels = "z score", 
              binwidth = 1, yLabel = "Number of region sets", 
              plotTitles = paste0("Distribution of region set z scores, ", colsToAnnotate),
              ggExpr = "+xlim(-3,40)")


topRSInd = rsRankingIndex(rsScores = rsZScores, signalCol = colsToAnnotate)


# get top region sets for each colsToAnnotate based on z score
topRSZAnnoList = list()
for (i in seq_along(colsToAnnotate)) {
    topRSZAnnoList[[i]] = data.frame(rsName=rsZScores$rsName[topRSInd[1:20, colsToAnnotate[i]]], 
                                     rsDescription=rsZScores$rsDescription[topRSInd[1:20, colsToAnnotate[i]]])
    names(topRSZAnnoList[[i]]) <- paste0(names(topRSZAnnoList[[i]]), "_", colsToAnnotate[i])
}

write.csv(topRSZAnnoList, file = paste0(sheetsDir, "topRSPermZScores", dataID, ".csv"))

# get top region sets based on p value








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
