.analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)
.plotSubdir = paste0(plotSubdir, "StatsPlots", .analysisID, "/")
if (!dir.exists(ffPlot(.plotSubdir))) {
    dir.create(ffPlot(.plotSubdir))
}

# load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", nPerm, "_", variationMetric,
#                    "_", dataID, ".RData")))
# load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", dataID, ".RData")))

# remove region sets that had no overlap
keepInd = apply(rsPermScores[[1]], MARGIN = 1, FUN = function(x) !any(is.na(x)))

nullDistList = lapply(X = 1:nrow(rsPermScores[[1]]),
                      FUN = function(x) extractNullDist(resultsList=rsPermScores, rsInd = x))

# screen out region sets with no overlap
nullDistList = nullDistList[keepInd]
realRSScores = realRSScores[keepInd, ]


# just an example of the null distributions for a single region set (arbitrarily rs1)
multiNiceHist(file = ffPlot(paste0(.plotSubdir, "nullDistRS1", .analysisID, ".pdf")), dataDF = nullDistList[[1]], 
              colsToPlot = colsToAnnotate, xLabels = "COCOA score (absolute correlation)", 
              binwidth = 0.001, yLabel = "Number of permutation results", 
              plotTitles = paste0("Null distribution of COCOA scores, ", colsToAnnotate),
              ggExpr = "+xlim(0,0.25)")


simpleCache(paste0("permPValsUncorrected", .analysisID), {
    rsPVals = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
                          calcCols=colsToAnnotate, whichMetric = "pval")
    rsPVals
}, recreate = recreate, reload = TRUE, assignToVariable = "rsPVals")
multiNiceHist(file = ffPlot(paste0(.plotSubdir, "pValDistUncorrected", .analysisID, ".pdf")), dataDF = rsPVals, 
              colsToPlot = colsToAnnotate, xLabels = "p-value", 
              binwidth = 0.005, yLabel = "Number of region sets", 
              plotTitles = paste0("Distribution of region set p-values, ", colsToAnnotate),
              ggExpr = paste0("+xlim(0,1)+ylim(0, ", nrow(rsPVals), ")"))

simpleCache(paste0("permZScores", .analysisID), { 
    rsZScores = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
                            calcCols=colsToAnnotate, whichMetric = "zscore")
    rsZScores
    
}, recreate = recreate, reload=TRUE, assignToVariable = "rsZScores")
# View(rsZScores[which(rsZScores$ > 7), ])

multiNiceHist(file = ffPlot(paste0(.plotSubdir, "zScoreDist", .analysisID, ".pdf")), dataDF = rsZScores, 
              colsToPlot = colsToAnnotate, xLabels = "z score", 
              binwidth = 1, yLabel = "Number of region sets", 
              plotTitles = paste0("Distribution of region set z scores, ", colsToAnnotate),
              ggExpr = "+xlim(-3,40)")


#topRSInd = rsRankingIndex(rsScores = rsZScores, signalCol = colsToAnnotate)

#################
# p-values based on fitted gamma distributions

# simpleCache(paste0("rsScore_", dataID), assignToVariable = "realRSScores")
correctionMethod = "BH" # input to p.adjust
gPValDF = getGammaPVal(scores = realRSScores[, colsToAnnotate, drop=FALSE], nullDistList = nullDistList, method = "mme", realScoreInDist = TRUE)
gPValDF = apply(X = gPValDF, MARGIN = 2, FUN = function(x) p.adjust(p = x, method = correctionMethod))
gPValDF = cbind(gPValDF, realRSScores[, colnames(realRSScores)[!(colnames(realRSScores) %in% colsToAnnotate)]])

simpleCache(paste0("permPValsCorrected", .analysisID), {
    gPValDF
}, recreate = recreate, reload = TRUE)

multiNiceHist(file = ffPlot(paste0(.plotSubdir, "corrected", correctionMethod, 
                                   "PValDist", .analysisID, ".pdf")), dataDF = gPValDF[colsToAnnotate],
              colsToPlot = colsToAnnotate, xLabels = "p-value",
              binwidth = .01, boundary = 0, yLabel = "Number of region sets",
              plotTitles = paste0("Distribution of region set p-values, ", colsToAnnotate),
              ggExpr = "+ylim(0, 2270) + xlim(0, 1)")
multiNiceHist(file = ffPlot(paste0(.plotSubdir, "corrected", correctionMethod, 
                                   "PValLog10Dist", .analysisID, ".pdf")), dataDF = -log10(gPValDF[colsToAnnotate]),
              colsToPlot = colsToAnnotate, xLabels = "p-value (-log10)",
              binwidth = 1, boundary = 0, yLabel = "Number of region sets",
              plotTitles = paste0("Distribution of region set p-values, ", colsToAnnotate),
              ggExpr = "+ylim(0, 2270)")

# p val cutoffs
sigCutoff = 0.05
trendCutoff = 0.1

# sort region sets according to p-value/z-score groups (significant, trending 
# toward significant, nonsignificant), then by average correlation score
pValGroups = as.data.frame(gPValDF)[, colsToAnnotate, drop=FALSE]

sigInd = pValGroups <= sigCutoff
notSigInd = pValGroups > trendCutoff
trendInd = (pValGroups <= trendCutoff) & (pValGroups > sigCutoff)

pValGroups[sigInd] = 1
pValGroups[notSigInd] = -1
pValGroups[trendInd] = 0


pValGroups = cbind(pValGroups, gPValDF[, c("rsName", "rsDescription")])
colnames(pValGroups) = paste0(colnames(pValGroups), "_PValGroup")

pRankedScores = cbind(realRSScores, pValGroups)
# View(dplyr::arrange(pRankedScores, desc(LF4_Z), desc(LF4)))
pRankedScores$index = 1:nrow(pRankedScores)
gPValDF2 = as.data.frame(gPValDF)[, colsToAnnotate, drop=FALSE]
colnames(gPValDF2) <- paste0(colnames(gPValDF2), "_PVal")
pRankedScores = cbind(pRankedScores, gPValDF2)


simpleCache(paste0("pRankedScores", .analysisID), {
  pRankedScores  
}, recreate = recreate)
# get top region sets for each colsToAnnotate based on p val
topRSZAnnoList = list()
topRSN = 50 # this many top RS for each colsToAnnotate
for (i in seq_along(colsToAnnotate)) {
    
    theseTopInd = dplyr::arrange(pRankedScores, 
                                 desc(get(paste0(colsToAnnotate[i], "_PValGroup"))), 
                                 desc(get(colsToAnnotate[i])))$index[1:topRSN]
    thesePValRanks = order(pRankedScores[, paste0(colsToAnnotate[i], "_PVal")], decreasing = FALSE)
    pRankedScores$index[thesePValRanks]
    topRSZAnnoList[[i]] = data.frame(pRankedScores[theseTopInd, c("rsName", "rsDescription", colsToAnnotate[i],
                                                                 paste0(colsToAnnotate[i], "_PValGroup"),  
                                                                 paste0(colsToAnnotate[i], "_PVal"),
                                                                 "signalCoverage", "regionSetCoverage", 
                                                                 "totalRegionNumber", "meanRegionSize")])
    
    names(topRSZAnnoList[[i]]) <- paste0(colsToAnnotate[i], "_", c("rsName", "rsDescription", "rsScore",
                                                                   "PValGroup", "pVal", "signalCoverage", "regionSetCoverage", 
                                                                   "totalRegionNumber", "meanRegionSize"))
}

write.csv(topRSZAnnoList, file = ffSheets(paste0("topRSPermpVals", .analysisID, ".csv")), row.names = FALSE)

#############################################################################
################
# get top region sets

# sort region sets according to p-value/z-score groups (significant, trending 
# toward significant, nonsignificant), then by average correlation score
zScoreGroups = as.data.frame(rsZScores)[, colsToAnnotate, drop=FALSE]
# 4.078 ~ abs(qnorm(0.05/2200))
sigCutoff = 4.078
# 1.645 ~ abs(qnorm(0.05))
trendCutoff = 1.645

zScoreGroups[zScoreGroups < trendCutoff] = -1
zScoreGroups[(zScoreGroups >= trendCutoff) & (zScoreGroups < sigCutoff)] = 0
zScoreGroups[zScoreGroups >= sigCutoff] = 1

zScoreGroups = cbind(zScoreGroups, rsZScores[, c("rsName", "rsDescription")])
colnames(zScoreGroups) = paste0(colnames(zScoreGroups), "_ZGroup")

zRankedScores = cbind(realRSScores, zScoreGroups)
# View(dplyr::arrange(realRSScores, desc(LF4_Z), desc(LF4)))
zRankedScores$index = 1:nrow(zRankedScores)
rsZScoresDF = as.data.frame(rsZScores)[, colsToAnnotate, drop=FALSE]
colnames(rsZScoresDF) <- paste0(colnames(rsZScoresDF), "_ZScore")
zRankedScores = cbind(zRankedScores, rsZScoresDF)

# simpleCache(paste0("zRankedScores", .analysisID), {
#     zRankedScores  
# })

# get top region sets for each colsToAnnotate based on z score
topRSZAnnoList = list()
topRSN = 50 # this many top RS for each colsToAnnotate
for (i in seq_along(colsToAnnotate)) {
    
    theseTopInd = dplyr::arrange(zRankedScores, 
                                 desc(get(paste0(colsToAnnotate[i], "_ZGroup"))), 
                                 desc(get(colsToAnnotate[i])))$index[1:topRSN]
    topRSZAnnoList[[i]] = data.frame(zRankedScores[theseTopInd, c("rsName", "rsDescription", colsToAnnotate[i],
                                                                 paste0(colsToAnnotate[i], "_ZGroup"),  
                                                                 paste0(colsToAnnotate[i], "_ZScore"),
                                                                 "signalCoverage", "regionSetCoverage", 
                                                                 "totalRegionNumber", "meanRegionSize")])
    
    names(topRSZAnnoList[[i]]) <- paste0(colsToAnnotate[i], "_", c("rsName", "rsDescription", "rsScore",
                                                                   "zGroup", "zScore", "signalCoverage", "regionSetCoverage", 
                                                                   "totalRegionNumber", "meanRegionSize"))
}

write.csv(topRSZAnnoList, file = ffSheets(paste0("topRSPermZScores", .analysisID, ".csv")), row.names = FALSE)
