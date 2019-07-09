# permutation test by shuffling sample labels

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "18-truePerm"


set.seed(1234)
nPerm = 500

######################################################################

# assigns methylMat, signalCoord, latentFactors to global environment
loadMOFAData()
genomicSignal = methylMat

# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg19")
simpleCache(paste0("rsScore_Cor_", "CLL196", "MOFA"), assignToVariable = "realRSScores")
row.names(realRSScores) = realRSScores$rsName
sharedRSNames = names(GRList)[names(GRList) %in% realRSScores$rsName]
GRList = GRList[sharedRSNames]
realRSScores = realRSScores[sharedRSNames, ]

#####################################################################
# functions for parallel permutations

# @param genomicSignal columns of dataMat should be samples/patients, rows should be genomic signal
# (each row corresponds to one genomic coordinate/range)
# @param sampleLabels Rows should be samples, columns should be "features" 
# (whatever you want to get correlation with: eg PC scores),
# all columns in featureMat will be used (subset when passing to function
# in order to not use all columns)
# @param calcCols character. the columns for which to calculate
# correlation and then to run COCOA on
corPerm <- function(randomInd, genomicSignal, 
                    signalCoord, GRList, calcCols,
                    sampleLabels) {
    
    # reorder the sample labels
    sampleLabels = sampleLabels[randomInd, ]
    
    # calculate correlation
    featureLabelCor = createCorFeatureMat(dataMat = genomicSignal, 
                                          featureMat = sampleLabels, 
                                          centerDataMat = TRUE, 
                                          centerFeatureMat = TRUE)
    
    # run COCOA
    thisPermRes = runCOCOA(loadingMat=featureLabelCor, 
                           signalCoord=signalCoord, GRList=GRList, 
                           PCsToAnnotate = calcCols, 
                           scoringMetric = "regionMean", verbose = TRUE)
    
    # return
    return(thisPermRes)
    
}
# This function will take a list of results of permutation tests that included
# many region sets and return a data.frame/data.table with the null
# distribution for a single region set (row)
# @param resultsList each item in the list is a data.frame, one item for
# each permutation with the results of that permutation. Each row in the 
# data.frame is a region set. Rows in all the data.frames should be
# in the same order.
# @param rsInd numeric. The row number for the region set of interest.
extractNullDist <- function(resultsList, rsInd) {
    rowList = lapply(resultsList, FUN = function(x) x[rsInd, ])
    rsNullDist = rbindlist(rowList)
    return(rsNullDist)
}


####################################################################

colsToAnnotate = paste0("LF", c(1:3, 5:7, 9))

indList = list()
# generate random indices for shuffling of samples
for (i in 1:nPerm) {
    indList[[i]] = sample(1:nrow(latentFactors), nrow(latentFactors))
}
# hist(sapply(indList, function(x) x[1]))


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
            genomicSignal=methylMat, 
            signalCoord=signalCoord, 
            GRList=GRList, 
            calcCols=colsToAnnotate,
            sampleLabels=latentFactors)
    message(i)
    if ((i %% 50) == 0) {
        save(rsPermScores, file = ffProc(paste0("COCOA_paper/RCache/rsPermScores_", scriptID, ".RData")))
    }
    
}

save(rsPermScores, file = ffProc(paste0("COCOA_paper/RCache/rsPermScores_", scriptID, ".RData")))

# get score null distribution for each region set
# one null distribution for each region set
nullDistList = lapply(X = seq_along(GRList),
                      FUN = function(x) extractNullDist(resultsList=rsPermScores, rsInd = x))
simpleCache("nullDistListMOFACor196", {
    nullDistList
})

hist(nullDistList[[1]]$LF1)

rsPVals = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
                      calcCols=colsToAnnotate, whichMetric = "pval")
hist(rsPVals$LF1, breaks = seq(0, 1, by=0.005))
hist(rsPVals$LF2, breaks = seq(0, 1, by=0.005))
hist(rsPVals$LF9, breaks = seq(0, 1, by=0.005))
hist(rsPVals$LF3)
hist(rsPVals$LF7)

View(rsPVals[which(rsPVals$LF3 < 0.01), ])

rsZScores = getPermStat(rsScores=realRSScores, nullDistList=nullDistList, 
                        calcCols=colsToAnnotate, whichMetric = "zscore")
View(rsZScores[which(rsZScores$LF7 > 25), ])
