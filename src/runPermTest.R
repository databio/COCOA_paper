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
#' @param variationMetric character. Either "cor" (correlation), "pcor" (partial
#' correlation), or "cov" (covariation)
#' 
# # for visualization
#' @param realRSScores


if (!exists("variationMetric")) {
    variationMetric = "cor"
}


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
                                sampleLabels=sampleLabels,
                                variationMetric = variationMetric)
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

#################
# simpleCache(paste0("rsScore_", dataID), assignToVariable = "realRSScores")
# p-values based on fitted gamma distributions
gPValDF = getGammaPVal(scores = realRSScores[, colsToAnnotate], nullDistList = nullDistList)
gPValDF = cbind(gPValDF, realRSScores[, colnames(realRSScores)[!(colnames(realRSScores) %in% colsToAnnotate)]])

multiNiceHist(file = ffPlot(paste0(plotSubdir, "pValDist", dataID, ".pdf")), dataDF = -log10(gPValDF[colsToAnnotate]),
              colsToPlot = colsToAnnotate, xLabels = "p-value",
              binwidth = 1, boundary = 0, yLabel = "Number of region sets",
              plotTitles = paste0("Distribution of region set p-values (-log10), ", colsToAnnotate),
              ggExpr = "+ylim(0, 2270)")

# p val cutoffs
sigCutoff = 0.05 / nrow(realRSScores)
trendCutoff = sigCutoff * 10

# sort region sets according to p-value/z-score groups (significant, trending 
# toward significant, nonsignificant), then by average correlation score
pValGroups = as.data.frame(gPValDF)[, colsToAnnotate]

sigInd = pValGroups <= sigCutoff
notSigInd = pValGroups > trendCutoff
trendInd = (pValGroups <= trendCutoff) & (pValGroups > sigCutoff)

pValGroups[sigInd] = 1
pValGroups[notSigInd] = -1
pValGroups[trendInd] = 0


pValGroups = cbind(pValGroups, gPValDF[, c("rsName", "rsDescription")])
colnames(pValGroups) = paste0(colnames(pValGroups), "_PValGroup")

realRSScores = cbind(realRSScores, pValGroups)
# View(dplyr::arrange(realRSScores, desc(LF4_Z), desc(LF4)))
realRSScores$index = 1:nrow(realRSScores)
gPValDF2 = as.data.frame(gPValDF)[, colsToAnnotate]
colnames(gPValDF2) <- paste0(colnames(gPValDF2), "_PVal")
realRSScores = cbind(realRSScores, gPValDF2)

# get top region sets for each colsToAnnotate based on p val
topRSZAnnoList = list()
topRSN = 50 # this many top RS for each colsToAnnotate
for (i in seq_along(colsToAnnotate)) {
    
    theseTopInd = dplyr::arrange(realRSScores, 
                                 desc(get(paste0(colsToAnnotate[i], "_PValGroup"))), 
                                 desc(get(colsToAnnotate[i])))$index[1:topRSN]
    thesePValRanks = order(realRSScores[, paste0(colsToAnnotate[i], "_PVal")], decreasing = FALSE)
    realRSScores$index[thesePValRanks]
    topRSZAnnoList[[i]] = data.frame(realRSScores[theseTopInd, c("rsName", "rsDescription", colsToAnnotate[i],
                                                                 paste0(colsToAnnotate[i], "_PValGroup"),  
                                                                 paste0(colsToAnnotate[i], "_PVal"),
                                                                 "signalCoverage", "regionSetCoverage", 
                                                                 "total_region_number", "mean_region_size")])
    
    names(topRSZAnnoList[[i]]) <- paste0(colsToAnnotate[i], "_", c("rsName", "rsDescription", "rsScore",
                                                                   "PValGroup", "pVal", "signalCoverage", "regionSetCoverage", 
                                                                   "total_region_number", "mean_region_size"))
}

write.csv(topRSZAnnoList, file = paste0(sheetsDir, "topRSPermpVals", dataID, ".csv"), row.names = FALSE)

################
# get top region sets

# sort region sets according to p-value/z-score groups (significant, trending 
# toward significant, nonsignificant), then by average correlation score
zScoreGroups = as.data.frame(rsZScores)[, colsToAnnotate]
# 4.078 ~ abs(qnorm(0.05/2200))
sigCutoff = 4.078
# 1.645 ~ abs(qnorm(0.05))
trendCutoff = 1.645

zScoreGroups[zScoreGroups < trendCutoff] = -1
zScoreGroups[(zScoreGroups >= trendCutoff) & (zScoreGroups < sigCutoff)] = 0
zScoreGroups[zScoreGroups >= sigCutoff] = 1

zScoreGroups = cbind(zScoreGroups, rsZScores[, c("rsName", "rsDescription")])
colnames(zScoreGroups) = paste0(colnames(zScoreGroups), "_ZGroup")

realRSScores = cbind(realRSScores, zScoreGroups)
# View(dplyr::arrange(realRSScores, desc(LF4_Z), desc(LF4)))
realRSScores$index = 1:nrow(realRSScores)
rsZScoresDF = as.data.frame(rsZScores)[, colsToAnnotate]
colnames(rsZScoresDF) <- paste0(colnames(rsZScoresDF), "_ZScore")
realRSScores = cbind(realRSScores, rsZScoresDF)

# get top region sets for each colsToAnnotate based on z score
topRSZAnnoList = list()
topRSN = 50 # this many top RS for each colsToAnnotate
for (i in seq_along(colsToAnnotate)) {
    
    theseTopInd = dplyr::arrange(realRSScores, 
                                 desc(get(paste0(colsToAnnotate[i], "_ZGroup"))), 
                                 desc(get(colsToAnnotate[i])))$index[1:topRSN]
    topRSZAnnoList[[i]] = data.frame(realRSScores[theseTopInd, c("rsName", "rsDescription", colsToAnnotate[i],
                                                                 paste0(colsToAnnotate[i], "_ZGroup"),  
                                                                 paste0(colsToAnnotate[i], "_ZScore"),
                                                                 "signalCoverage", "regionSetCoverage", 
                                                                 "total_region_number", "mean_region_size")])
        
    names(topRSZAnnoList[[i]]) <- paste0(colsToAnnotate[i], "_", c("rsName", "rsDescription", "rsScore",
                                           "zGroup", "zScore", "signalCoverage", "regionSetCoverage", 
                                           "total_region_number", "mean_region_size"))
}

write.csv(topRSZAnnoList, file = paste0(sheetsDir, "topRSPermZScores", dataID, ".csv"), row.names = FALSE)


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
