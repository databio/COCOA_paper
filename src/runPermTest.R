# permutation test by shuffling sample labels
# expected inputs: 
#' permutation test by shuffling sample labels
#' 
#' For reproducibility, set seed with 'set.seed()' function before running.
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
#' @param useSimpleCache logical. 
#' @param cacheDir character.
#' @param correctionMethod character. P value correction method. For acceptable 
#' arguments see ?p.adjust() (method parameter) 
#' @param resultType character. "pval" or "zscore"
#' @param ... character. Optional additional arguments for simpleCache
#' 
# # for visualization
#' @param realRSScores
#' @return Returns a list where each item is a data.frame of COCOA results 
#' from a separate permutation

#' runCOCOAPerm <- function(genomicSignal, 
#'                                  signalCoord, 
#'                                  GRList, 
#'                                  signalCol=c("PC1", "PC2"),
#'                                  signalCoordType = "default", 
#'                                  scoringMetric="default",
#'                                  variationMetric="cor",
#'                                  nPerm=300,
#'                                  useSimpleCache=TRUE,
#'                                  cacheDir=getwd(),
#'                                  dataID="tmp", ...) {
#'                                  
#' 
# indList <- list()
# # generate random indices for shuffling of samples
# for (i in 1:nPerm) {
#     indList[[i]] <- sample(1:nrow(sampleLabels), nrow(sampleLabels))
# }
# 
# # # plug random indices into parallelized function
# # rsPermScores <- COCOA:::lapplyAlias(X= indList, FUN=function(x) corPerm(randomInd=x, 
# #                                                         genomicSignal=methylMat, 
# #                                                         signalCoord=signalCoord, 
# #                                                         GRList=GRList, 
# #                                                         calcCols=colsToAnnotate,
# #                                                         sampleLabels=latentFactors))
# 
# # create the main permutation cache
# simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), {
#     
#     rsPermScores <- list()
#     for (i in seq_along(indList)) {
#         # for (i in (length(rsPermScores) + 1):nPerm) {
#         onePermCacheName <- paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID, "_Cache", i)
#         # create sub caches, one for each permutation
#         simpleCache(onePermCacheName, cacheSubDir = paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), {
#             
#             tmp <- corPerm(randomInd=indList[[i]], 
#                            genomicSignal=genomicSignal, 
#                            signalCoord=signalCoord, 
#                            GRList=GRList, 
#                            calcCols=colsToAnnotate,
#                            sampleLabels=sampleLabels,
#                            variationMetric = variationMetric)
#             message(i) # must be ahead of object that is saved as cache, not after
#             tmp
#             
#         })
#     }
#     
#     # combining all individual permutations/caches into one object
#     for (i in seq_along(indList)) {
#         onePermCacheName <- paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID, "_Cache", i)
#         rsPermScores[[i]] <- get(onePermCacheName)
#     }
#     rsPermScores
#     
# }, assignToVariable="rsPermScores")
#' 
#' .analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)
# .plotSubdir = paste0(plotSubdir, "StatsPlots", .analysisID, "/")
# if (!dir.exists(ffPlot(.plotSubdir))) {
#     dir.create(ffPlot(.plotSubdir))
# }
# 
# # remove region sets that had no overlap
# keepInd = apply(rsPermScores[[1]], MARGIN = 1, FUN = function(x) !any(is.na(x)))
# 
# nullDistList = lapply(X = 1:nrow(rsPermScores[[1]]),
#                       FUN = function(x) extractNullDist(resultsList=rsPermScores, rsInd = x))
# 
# # screen out region sets with no overlap
# nullDistList = nullDistList[keepInd]
# realRSScores = realRSScores[keepInd, ]
# 
# simpleCache(paste0("permPValsUncorrected", .analysisID), {
#     rsPVals = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
#                           calcCols=colsToAnnotate, whichMetric = "pval")
#     rsPVals
# }, recreate = TRUE, reload = TRUE)
# 
# simpleCache(paste0("permZScores", .analysisID), { 
#     rsZScores = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
#                             calcCols=colsToAnnotate, whichMetric = "zscore")
#     rsZScores
#     
# }, recreate = TRUE, reload=TRUE)
# 
# 
# #topRSInd = rsRankingIndex(rsScores = rsZScores, signalCol = colsToAnnotate)
# 
# #################
# # simpleCache(paste0("rsScore_", dataID), assignToVariable = "realRSScores")
# # p-values based on fitted gamma distributions
# correctionMethod = "BH" # input to p.adjust
# gPValDF = getGammaPVal(scores = realRSScores[, colsToAnnotate, drop=FALSE], nullDistList = nullDistList, method = "mme", realScoreInDist = TRUE)
# gPValDF = apply(X = gPValDF, MARGIN = 2, FUN = function(x) p.adjust(p = x, method = correctionMethod))
# gPValDF = cbind(gPValDF, realRSScores[, colnames(realRSScores)[!(colnames(realRSScores) %in% colsToAnnotate)]])
# 
# simpleCache(paste0("permPValsCorrected", .analysisID), {
#     gPValDF
# }, recreate = TRUE, reload = TRUE)
#' 
#' }
###############################################################################

if (!exists("variationMetric")) {
    variationMetric <- "cor"
}


indList <- list()
# generate random indices for shuffling of samples
for (i in 1:nPerm) {
    indList[[i]] <- sample(1:nrow(sampleLabels), nrow(sampleLabels))
}

# # plug random indices into parallelized function
# rsPermScores <- COCOA:::lapplyAlias(X= indList, FUN=function(x) corPerm(randomInd=x, 
#                                                         genomicSignal=methylMat, 
#                                                         signalCoord=signalCoord, 
#                                                         GRList=GRList, 
#                                                         calcCols=colsToAnnotate,
#                                                         sampleLabels=latentFactors))

# create the main permutation cache
simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), {
    
    rsPermScores <- list()
    for (i in seq_along(indList)) {
        # for (i in (length(rsPermScores) + 1):nPerm) {
        onePermCacheName <- paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID, "_Cache", i)
        # create sub caches, one for each permutation
        simpleCache(onePermCacheName, cacheSubDir = paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), {
            
            tmp <- corPerm(randomInd=indList[[i]], 
                          genomicSignal=genomicSignal, 
                          signalCoord=signalCoord, 
                          GRList=GRList, 
                          calcCols=colsToAnnotate,
                          sampleLabels=sampleLabels,
                          variationMetric = variationMetric)
            message(i) # must be ahead of object that is saved as cache, not after
            tmp
            
        })
    }
    
    # combining all individual permutations/caches into one object
    for (i in seq_along(indList)) {
        onePermCacheName <- paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID, "_Cache", i)
        rsPermScores[[i]] <- get(onePermCacheName)
    }
    rsPermScores
    
}, assignToVariable="rsPermScores")


source(ffProjCode("processPermResults.R"))

## renaming objects that were created manually to make them compatible with simpleCache
# rm(rsPermScores)
# rm(ret)
# load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", nPerm, "_", variationMetric, "_", dataID, ".RData")))
# ret <- rsPermScores
# save(ret, file =ffProc(paste0("COCOA_paper/RCache/rsPermScores_", nPerm, "_", variationMetric, "_", dataID, ".RData")))


