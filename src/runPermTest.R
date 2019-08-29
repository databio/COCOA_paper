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

simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), {
    rsPermScores = list()
    #for (i in seq_along(indList)) {
    for (i in (length(rsPermScores) + 1):nPerm) {
        
        rsPermScores[[i]] = corPerm(randomInd=indList[[i]], 
                                    genomicSignal=genomicSignal, 
                                    signalCoord=signalCoord, 
                                    GRList=GRList, 
                                    calcCols=colsToAnnotate,
                                    sampleLabels=sampleLabels,
                                    variationMetric = variationMetric)
        message(i)
        save(rsPermScores, file = ffProc(paste0("COCOA_paper/RCache/rsPermScores_", 
                                                dataID, ".RData")))
        # if ((i %% 50) == 0) {
        #     save(rsPermScores, file = ffProc(paste0("COCOA_paper/RCache/rsPermScores_", 
        #                                             dataID, ".RData")))
        # }
        
    }
    
    save(rsPermScores, file = ffProc(paste0("COCOA_paper/RCache/rsPermScores_", 
                                            variationMetric, "_", 
                                            dataID, ".RData")))
    rsPermScores
    
}, assignToVariable="rsPermScores")


source(ffProjCode("processPermResults.R"))

## renaming objects that were created manually to make them compatible with simpleCache
# rm(rsPermScores)
# rm(ret)
# load(ffProc(paste0("COCOA_paper/RCache/rsPermScores_", nPerm, "_", variationMetric, "_", dataID, ".RData")))
# ret = rsPermScores
# save(ret, file =ffProc(paste0("COCOA_paper/RCache/rsPermScores_", nPerm, "_", variationMetric, "_", dataID, ".RData")))
