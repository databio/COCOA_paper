# test/visualize different scoring metrics and normalizers
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
scriptID = "00_testScoring"

############################# testing out various scoring methods
# not unit tests, just exploratory analysis

# testing out scoring metric
totalCpGs = 350000
cytosine_coverage = 1:(floor(totalCpGs/2))
sizeNormalizer = 1 / sqrt((1 / cytosine_coverage) - (1 / (totalCpGs - cytosine_coverage)))
plot(sizeNormalizer)
test = sizeNormalizer[2:(length(sizeNormalizer))] - sizeNormalizer[1:(length(sizeNormalizer)-1)]
plot(test)
plot(test[50000:100000])

# rank sum test
n = 1:100
mDist = lapply(n, function(x) rnorm(350000 - x))
wTest = mapply(FUN = function(x, y) wilcox.test(x = rep(0.67, x), y = y), x = n, y=mDist)
length(wTest)
# plot the test statistic
plot(unlist(wTest[(n * 7) - 6]))
# plot the p value
plot(unlist(wTest[(n * 7) - 4])[ 40:100])

n = 1:10
mDist = lapply(n, function(x) rnorm(1600000 - x)) # 1.6 million CpGs
mDist2 = lapply(n, function(x) rnorm(x, mean = 0.5))
wTest = mapply(FUN = function(x, y) wilcox.test(x = x, y = y), x = mDist2, y=mDist)
# plot the test statistic
plot(unlist(wTest[(n * 7) - 6]))
# plot the p value
plot(unlist(wTest[(n * 7) - 4])[ 40:100])

# combining raw score and test statistic approach
rsEn_old[, rsIndex := 1:nrow(rsEn_old)]
rsEn_new[, rsIndex := 1:nrow(rsEn_new)]
oldEn = rsEn_old
oldEn = oldEn[order(PC1, decreasing = TRUE)]
oldEn[, PC1 := 1:nrow(oldEn)]
oldEn = oldEn[order(rsIndex), ]

newEn = rsEn_new
newEn = newEn[order(PC1, decreasing = TRUE)]
newEn[, PC1 := 1:nrow(newEn)]
newEn = newEn[order(rsIndex), ]

newEn[, newOrder := PC1 + oldEn$PC1]
View(newEn[order(newOrder), ])

# pooled standard deviation
SD1 = 50
SD2 = 1:100
poolSD = sqrt((SD1^2 + SD2^2)/ 2)
plot(poolSD)


regionGR = GRList[[as.numeric(rsEnSortedInd[1,1])]]

# magnitude of a single vector (distance from zero)
vecDist = function(dataVec) {
    zVec = rep(0, length(dataVec))
    return(as.numeric(dist(rbind(zVec, dataVec))))
}

# get angle between PC and only loading vals for CpGs in a given region set
subsetPCAngle <- function(regionGR, loadingMat, coordinateDT, PCofInterest=paste0("PC", 1:3)) {
    #  subsetInd = 
    # newColNames <- paste0(PCofInterest, "_subset")
    
    coordGR = MIRA:::dtToGr(coordinateDT)
    olList = findOverlaps(query = regionGR, subject = coordGR)
    # regionHitInd = sort(unique(queryHits(olList)))
    cytosineHitInd = sort(unique(subjectHits(olList)))
    subsetLoad = loadingMat[, PCofInterest]
    subsetLoad[-cytosineHitInd, ] = 0
    # angle between two vector = (A * B) / ()
    subMag = apply(subsetLoad, 2, FUN = vecDist)
    dotProd = mapply(FUN = function(x, y) loadingMat[, x] * subsetLoad[, y], x = seq_along(PCofInterest), y=seq_along(PCofInterest))
    # sqrt(apply(dotProd, 2, sum)), could just take sqrt of subMag for dotProd
    dotProd = apply(dotProd, 2, sum)
    subAngle = acos(dotProd/ (apply(loadingMat[, PCofInterest], 2, vecDist) * subMag))
    # answer is same as taking acos(subMag)
    degAngle = subAngle * 180 / pi
}

subsetDist = function(dataVec, subsetInd) {
    
    pairMat = rbind(dataVec, dataVec)
    # set everything not in the subset to 0
    pairMat[2, -subsetInd] = 0
    return(as.numeric(dist(pairMat)))
}



nC = nrow(loadingMat)
lengthMat = t(cbind(loadingMat[, "PC1"], rep(0, nC)))
dist(lengthMat)

############################################################################
# testing wilcox scoring with confidence intervals

###########################################################
# reading in the region sets
# load LOLA database
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))

#################################################################
PCsToAnnotate = paste0("PC", 1:4)
rsScoreCacheName = paste0(scriptID, "_WilcoxConfInt657")
overwriteRSScoreResultsCaches = TRUE

simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "bigSharedC", reload = TRUE)

# GRList # from load_process_regions pipeline
signalCoord = bigSharedC$coordinates
allMPCAString = "allMPCA_657"
simpleCache(allMPCAString, assignToVariable = "mPCA", reload = TRUE)
loadingMat = mPCA$rotation
# use rsEnString to specify?
simpleCache("rsEnrichmentRSMean_657", assignToVariable = "rsScores", reload = TRUE)
# TODO make sure GRList and rsScores are both in the same order/with same data
names(GRList) <- paste0(rsScores$rsName, " : ", rsScores$rsDescription)
GRList = GRList[!is.na(rsScores$PC1)]
rsScores=rsScores[!is.na(rsScores$PC1), ]
rsName = rsName[!is.na(rsScores$PC1)]
rsDescription = rsDescription[!is.na(rsScores$PC1)]

# the region set database has more region sets than rsScores so indexing by
# name instead of row number
names(rsName) = rsName
names(GRList) = rsName
names(rsDescription) = rsName

rsEnSortedInd= rsRankingIndex(rsScores = rsScores, PCsToAnnotate = PCsToAnnotate)

rsInd = c(unique(as.numeric(as.matrix(rsEnSortedInd[1:10, PCsToAnnotate]))),
          unique(as.numeric(as.matrix(rsEnSortedInd[1001:1005, PCsToAnnotate]))),
          unique(as.numeric(as.matrix(rsEnSortedInd[(nrow(rsEnSortedInd)-4):nrow(rsEnSortedInd), PCsToAnnotate]))))
rsIndNames = rsScores$rsName[rsInd]
subGRList = GRList[rsIndNames]
subRSNames = rsName[rsIndNames]
subRSDescription = rsDescription[rsIndNames]


scoringMetric = "rankSum"

simpleCache(rsScoreCacheName, {
    rsScores = runCOCOA(loadingMat=loadingMat, 
                       signalCoord = signalCoord, 
                       GRList = subGRList, 
                       PCsToAnnotate = PCsToAnnotate, 
                       scoringMetric=scoringMetric,
                       verbose = TRUE, 
                       wilcox.conf.int = TRUE)
    rsScores$rsName = subRSNames
    rsScores$rsDescription= subRSDescription
    rsScores
}, recreate=overwriteRSScoreResultsCaches)

View(rsScore[order(rsScore$surv_at_daysCutoff, decreasing = TRUE), ])
hist(rsScore$surv_at_daysCutoff)

write.csv(x = rsScore, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/rsScore_brcaMethylPCR_", dataID, ".csv"),
          quote = FALSE, row.names = FALSE)




    




