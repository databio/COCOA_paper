# permutation test by shuffling sample labels

nCores = detectCores()
setOption("mc.cores"=nCores)

set.seed(1234)
nPerm = 16

######################################################################

# assigns methylMat, signalCoord, latentFactors to global environment
loadMOFAData()

# loads database of region sets 
# (assigns GRList, rsName, rsDescription to global environment)
loadGRList(genomeV="hg38")

#####################################################################


indList = list()
# generate random indices for shuffling of samples
for (i in 1:nPerm) {
    indList[[i]] = sample(1:5, 5)
}
# hist(sapply(indList, function(x) x[1]))


# plug random indices into parallelized function


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
colsToAnnotate = paste0("LF", c(1:3, 5:7, 9))

rsPermScores = COCOA:::lapplyAlias(X= indList, FUN=function(x) corPerm(randomInd=x, 
                                                        genomicSignal=methylMat, 
                                                        signalCoord=signalCoord, 
                                                        GRList=GRList, 
                                                        calcCols=colsToAnnotate,
                                                        sampleLabels=latentFactors))


# get score null distribution for each region set

extractNullDist <- function(resultsList, rsInd) {
    rowList = lapply(resultsList, FUN = function(x) x[rsInd, ])
    rsNullDist = rbindlist(rowList)
    return(rsNullDist)
}

# one null distribution for each region set
nullDistList = lapply(X = seq_along(GRList),
                      FUN = function(x) extractNullDist(resultsList=rsPermScores, rsInd = x))

