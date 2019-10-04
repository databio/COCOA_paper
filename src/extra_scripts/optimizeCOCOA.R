# make COCOA faster/more efficient

library(rbenchmark)
devtools::load_all(paste0(Sys.getenv("CODE"), "COCOA"))

# test how long COCOA takes at baseline

data("brcaMCoord1")
data("brcaLoadings1")
data("esr1_chr1")
data("nrf1_chr1")
data("brcaMethylData1")
data("brcaPCScores657")
sampleLabels = brcaPCScores657[colnames(brcaMethylData1), ]
sampleLabels$ER_Status = scale(as.numeric(sampleLabels$ER_Status),
                               center=TRUE, scale=FALSE)
signalCoord = brcaMCoord1
calcCols = "ER_Status"
GRList = GRangesList(esr1_chr1, nrf1_chr1)
# shuffling sample labels
randomInd = sample(1:nrow(sampleLabels), nrow(sampleLabels))

benchmark(corPerm(randomInd=randomInd, genomicSignal=brcaMethylData1,
        signalCoord=brcaMCoord1, GRList=GRangesList(esr1_chr1, nrf1_chr1),
        calcCols="ER_Status", sampleLabels=sampleLabels,
        variationMetric="cor"), replications = 50)


# replications elapsed relative user.self sys.self user.child sys.child
# 50  55.596        1    64.796        0          0         0

### after switching nested apply() for cor(matrix, matrix) call in createCorFeatureMat 
# replications elapsed relative user.self sys.self user.child sys.child
# 50  15.438        1    24.336    0.016          0         0

### after removing unnecessary copy: dataMat = data.table::copy(as.data.frame(t(dataMat)))
# replications elapsed relative user.self sys.self user.child sys.child
# 50  14.441        1    23.256    0.012          0         0

######################### optimizing runCOCOA ##################################
variationMetric = "cor"
absVal = TRUE
verbose= TRUE

featureLabelCor = createCorFeatureMat(dataMat = brcaMethylData1, 
                                      featureMat = sampleLabels[, "ER_Status", drop=FALSE], 
                                      centerDataMat = TRUE, 
                                      centerFeatureMat = TRUE,
                                      testType = variationMetric)[, "ER_Status", drop=FALSE]

benchmark(runCOCOA(signal=featureLabelCor, 
         signalCoord=signalCoord, GRList=GRList, 
         signalCol = calcCols, 
         scoringMetric = "default", verbose = verbose,
         absVal = absVal), replications = 50) 
# replications elapsed relative user.self sys.self user.child sys.child
# 50   3.162        1    12.324    0.004          0         0

########################## optimizing createCorFeatureMat ############################


benchmark(createCorFeatureMat(dataMat = brcaMethylData1, 
                        featureMat = sampleLabels, 
                        centerDataMat = TRUE, 
                        centerFeatureMat = TRUE,
                        testType = variationMetric), replications=50)

# replications elapsed relative user.self sys.self user.child sys.child
# 50  10.961        1    10.932    0.008          0         0

##### with centering = TRUE
# replications elapsed relative user.self sys.self user.child sys.child
# 50  11.237        1    11.224    0.012          0         0



############## other stuff in corPerm ################################

otherStuff <- function(randomInd, genomicSignal, 
                               signalCoord, GRList, calcCols,
                               sampleLabels, variationMetric = "cor", 
                               scoringMetric="default", verbose=TRUE,
                               absVal=TRUE) {
    # if vector is given, return error
    if (is.null(dim(sampleLabels))) {
        stop("`sampleLabels` should be a matrix or data.frame")
    }
    
    if (any(!(calcCols %in% colnames(sampleLabels)))) {
        stop("Not all specified columns are present in `sampleLabels`")
    }
    
    # subset to only calcCols
    sampleLabels = sampleLabels[, calcCols, drop=FALSE]
    
    # because names are dropped for a single column data.frame when indexing
    # single col data.frame is automatically converted to numeric
    featureNames = colnames(sampleLabels)
    # reorder the sample labels
    sampleLabels = data.frame(sampleLabels[randomInd, ])
    colnames(sampleLabels) = featureNames
}

benchmark(otherStuff(randomInd=randomInd, genomicSignal=brcaMethylData1,
                  signalCoord=brcaMCoord1, GRList=GRList,
                  calcCols="ER_Status", sampleLabels=sampleLabels,
                  variationMetric="cor"), replications = 50)

# replications elapsed relative user.self sys.self user.child sys.child
# 50   0.018        1      0.02        0          0         0
