## generate simulated DNA methylation data to compare COCOA and LOLA

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
library(curatedTCGAData)
library(TCGAutils)
library(bumphunter)
library(LOLA)
loadCOCOA()

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
plotSubdir = "50-Sim/"

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

plotWidth = 80
plotHeight = 80
plotUnits = "mm"

genomeV = "hg19"
dataID = "KIRC_simulated"

set.seed(1234)


##################################
# start from single healthy sample
loadGRList(genomeV = genomeV)

simpleCache(paste0("healthy_KIRC_methyl_", dataID) , {
    loadTCGAMethylation(cancerID = "KIRC")
    methylMat = methylList$methylProp
    signalCoord = methylList$coordinates
    
    sampleType = substr(colnames(methylMat), start = 14, stop = 15)
    # 01 is primary solid tumor, 11 is solid normal tissue, 05 is new primary tumor
    # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
    # https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
    normalSampleInd = (sampleType == "11")
    tumorSampleInd = (sampleType == "01") # exclude the extra sample 05
    methylMat = methylMat[, normalSampleInd]
    # now I only need patient ID
    colnames(methylMat) = substr(colnames(methylMat), start = 1, stop = 12)
    
    # average of all normal samples (160 samples)
    healthy = apply(X = methylMat, 1, mean)
    methylList$methylProp = healthy
    methylList
}, assignToVariable = "methylList")

signalCoord = methylList$coordinates
signalCoord = COCOA:::dtToGr(signalCoord)

simpleCache(paste0("both_", dataID, "_no_noise"), {
    healthy = methylList$methylProp
    
    # fake tumor profile, only changing DNA methylation in one region set
    tumor = healthy
    
    # find CpGs that overlap region set
    myGR = GRList[[c("Human_MDA-MB-231-Cells_ESR1,-DBDmut_E2-45min_Katzenellenbogen.bed")]] # 356 
    # # H1hesc_WE.bed
    fo = findOverlaps(query = myGR, subject = signalCoord)
    olCGInd = subjectHits(fo)
    tumor[olCGInd] = 0.1
    # set healthy sample to opposite of tumor in this region set
    healthy[olCGInd] = 0
    
    # make samples with varying proportion healthy vs "tumor"
    # name says what percent tumor
    percTumor = seq(10, 100, 10)
    percHealth = 100 - percTumor
    mixed = ((tumor %*% t(percTumor)) + (healthy %*% t(percHealth))) / 100
    colnames(mixed) = paste0("tumor", percTumor)
    
    both = matrix(data = healthy, nrow = length(healthy), ncol=10)
    colnames(both) = paste0("healthy", 1:10)
    # 10 healthy, 10 tumor
    both = cbind(both, mixed) 
}, assignToVariable = "both")



###############


# Function gets DMRs with bumphunter then runs LOLA on those
# @param signalCoord GRanges.
# @param GRList GrangesList for LOLA
# @param rsAnno data.frame. One row for each region set in GRList.
# @param targetVar numeric. Can be discrete or continuous. Used with bumphunter. 
dmrLOLA = function(genomicSignal, signalCoord, GRList, rsAnno, targetVar, dataID="") {
    
    both = genomicSignal
    sampleStatus = targetVar
    
    designMat = model.matrix(object = ~ sampleStatus)

    
    simpleCache(paste0("tumorDMR_", dataID), {
        tumorDMR =bumphunter(object = both, design=designMat, chr=COCOA:::grToDt(signalCoord)$chr, 
                             pos=COCOA:::grToDt(signalCoord)$start, coef=2, cutoff = 0.1, B=1000)# , type="Beta")
        tumorDMR
    }, assignToVariable = "tumorDMR")
    # View(tumorDMR$table)
    
    # LOLA
    if (nrow(tumorDMR$table[tumorDMR$table$fwer <=0.05, c("chr", "start", "end")]) == 0) {
        return(NULL)
    }
    dmrGR = COCOA:::dtToGr(tumorDMR$table[tumorDMR$table$fwer <=0.05, c("chr", "start", "end")])
    uniGR = resize(x = signalCoord, width = 2000, fix = "center")
    uniGR = reduce(uniGR)
    names(GRList) = NULL
    regionSetDB = list()
    regionSetDB$regionAnno = rsAnno
    regionSetDB$regionGRL = GRList
    lResults = LOLA::runLOLA(userSets = dmrGR, userUniverse = uniGR, regionDB = regionSetDB)
    # View(arrange(lResults, desc(oddsRatio)))
    # View(arrange(lResults, desc(pValueLog)))
    # sum(lResults$pValueLog < 2) / nrow(lResults)  
    
    return(lResults)
}

##################################
# compare p-value approximation to actual p-value

###########
simpleCache("mixedGRList", {
    
    fo = findOverlaps(query = myGR, subject = signalCoord)
    olCGInd = subjectHits(fo)
    # create test region sets
    rsOlInd = queryHits(fo)
    myGR = myGR[unique(rsOlInd)]
    rs100 = myGR
    groupSize = floor(length(unique(rsOlInd))/ 20) 
    realRSSize = length(myGR)
    
    # make a region set of random regions (that do not overlap with GR of interest)
    randomGR = resize(signalCoord[-olCGInd][1:realRSSize * 100], width= 500, fix="center")
    mixedGRList = list()
    
    for (i in 1:20) {
        mixedGRList[[i]] = c(myGR[(1+(i-1)*groupSize):length(myGR)], 
                             randomGR[0:((i-1)*groupSize)])
    }
    # only include a few real regions and almost all random regions
    for (j in 1:5) {
        mixedGRList[[i + j]] = c(myGR[(length(myGR)-5+j):length(myGR)], 
                                 randomGR[0:(length(myGR)-6+j)])
    }
    mixedGRList[[i+j+1]] = randomGR
    
    
    mixedGRList = GRangesList(mixedGRList)
    mixNames = c(paste0("gr", seq(100, 5, -5)), paste0("gr", seq(5,0, -1), "RealRegions"))
    names(mixedGRList) = mixNames
    mixedGRList
    
}, assignToVariable = "mixedGRList")
mixNames = names(mixedGRList)
##################################
# running simulated data with noise, compare COCOA to LOLA
set.seed(1234)

noise = matrix(rnorm(n = nrow(both) * ncol(both), mean = 0, sd = 0.05), 
               nrow = nrow(both), ncol = ncol(both))
bothHigh = both +  noise
bothHigh[bothHigh < 0] = 0
bothHigh[bothHigh > 1] = 1
bothHPCA = prcomp(t(bothHigh))$x
plot(bothHPCA[, c(1, 2)])
wApplyFun = function(x, pcScores) {
    wilcox.test(pcScores[paste0("healthy", 1:10), x], pcScores[paste0("tumor", 1:10 *10), x])$p.value
}
wRes = sapply(X = paste0("PC", 1:10), FUN = function(x) wApplyFun(x,bothHPCA))
wRes

noise = matrix(rnorm(n = nrow(both) * ncol(both), mean = 0, sd = 0.025), 
               nrow = nrow(both), ncol = ncol(both))
bothLow = both +  noise
bothLow[bothLow < 0] = 0
bothLow[bothLow > 1] = 1
bothLPCA = prcomp(t(bothLow))$x
wRes = sapply(X = paste0("PC", 1:10), FUN = function(x) wApplyFun(x, bothLPCA))
wRes

########### figure with samples ordered by PC score
# order samples by PC score, color by ER status
pcScores = bothHPCA
erStatusDF = cbind(apply(X = pcScores[, paste0("PC", 1:2)], 
                         MARGIN = 2, FUN = function(x) order(order(x, decreasing = FALSE))), 
                   as.data.frame(patientMetadata[row.names(pcScores) , ])[, c("subject_ID", "ER_status")])
erStatusDF = pivot_longer(erStatusDF, cols = c("PC1", "PC2", "PC3", "PC4"), names_to = "PC", values_to = "rank")    

erStatusDF$barHeight = rep(1, nrow(erStatusDF))
erStatusDF = arrange(erStatusDF, PC, rank) 

erStatusPlot = ggplot(data = erStatusDF, mapping = aes(x=rank, y=barHeight, group=PC)) + 
    geom_col(aes(col=ER_status)) + scale_color_discrete(breaks=c("Positive", "Negative")) +
    xlab("Samples ordered by PC score") + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                                axis.title.y = element_blank(), axis.line = element_blank())
erStatusPlot

ggplot2::ggsave(filename=ffPlot(paste0(plotSubdir,"/orderedERStatus.svg")), 
                plot = erStatusPlot, device = "svg", height = plotHeight / 2, 
                width = plotWidth, units = plotUnits)




t######################################################
# targetVar = bothLPCA[, paste0("PC", 1:10)]
# 
# nreps = 10000
# m = replicate(nreps, sample(seq_len(nrow(targetVar))))
# 
# targetVar
# targetVar[m[,1],]
# 
# mt = matrix(targetVar[m], nrow=nrow(targetVar))
# mt
# dim(mt)
# colnames(mt) = paste0("rep", seq_len(nreps))
# 
# dim(bothLow)
# 
# res = cor(t(brcaMethylData1), mt)
# dim(res)
# res
# 
#############

# a = findOverlaps(query = myGR, mixedGRList[[10]])
# 


######################################################

###################################
# compare empirical p-value to gamma pval after adding noise

# filter data to only relevant CpGs and coordinates to speed things up
resList = getOLRegions(GRList = mixedGRList, intGR = signalCoord)
signalCoordTmp = signalCoord[resList$intGROverlapInd]
bothHighTmp = bothHigh[resList$intGROverlapInd, ]

subCache = paste0(getCacheDir(), "/gammaSimulations")
signalCol = c("PC1", "PC2")
dir.create(subCache)
smallChange = runCOCOA(genomicSignal = bothHighTmp, signalCoord = signalCoordTmp, GRList = mixedGRList, 
                       signalCol = c("PC1", "PC2"), targetVar = bothHPCA, 
                       variationMetric = "cov", scoringMetric = "regionMean", 
                       absVal = TRUE, centerGenomicSignal = TRUE, centerTargetVar = TRUE)
smallChange = cbind(smallChange, rsName=mixNames)

# rm(list = c("methylList", "regionSetDB"))
# gc()
set.seed(1234)
setLapplyAlias(cores = 6)
a=runCOCOAPerm(nPerm = 100000, rsScores = smallChange[, signalCol], useSimpleCache = TRUE, 
               cacheDir = subCache, dataID = "simWithNoise05",
               genomicSignal = bothHighTmp, signalCoord = signalCoordTmp, GRList = mixedGRList, 
               signalCol = signalCol, targetVar = bothHPCA[, signalCol], 
               variationMetric = "cov", scoringMetric = "regionMean", 
               absVal = TRUE, centerGenomicSignal = TRUE, centerTargetVar = TRUE, recreate=FALSE)
View(a$empiricalPVals)
View(a$gammaPVal)
nullDistList = convertToFromNullDist(a$permRSScores)
hist(nullDistList[[1]]$PC1)
loadPermCaches = function(n, cacheDir, cacheString, assignToVariable) {
    rsPermScores = list()
    for (i in n) {
        simpleCache(paste0(cacheString, i), cacheDir = cacheDir, assignToVariable = "rsPermScoresTmp")
        rsPermScores[[i]] <- rsPermScoresTmp
    }
    assign(x =assignToVariable, value = rsPermScores, envir = parent.frame(n=1))
}
# loadPermCaches(n=1:100000, paste0(getCacheDir(), "/gammaSimulations/rsPermScores_20000Perm_cov_simWithNoise05/"), 
#                cacheString = "rsPermScores_20000Perm_cov_simWithNoise05_Cache", 
#                assignToVariable = "rsPermScores")
# nullDist = COCOA::convertToFromNullDist(rsPermScores)
# hist(nullDist[[26]]$PC1)
# rsPVals <- COCOA:::getPermStat(rsScores=smallChange, nullDistList=nullDist,
#                        signalCol=signalCol, whichMetric = "pval",
#                        testType="greater")
# View(rsPVals)
# length(rsPVals)
# hist(nullDist[[26]]$PC1)

# for (i in 1:100000) {
    # simpleCache(paste0("rsPermScores_1e+05Perm_cov_simWithNoise05_Cache", i),
    #             cacheDir = paste0(subCache, "/", "rsPermScores_1e+05Perm_cov_simWithNoise05/"),
    #             assignToVariable = "tmp")
#     tmp = tmp[, c("PC1", "PC2")]
#     simpleCache(paste0("rsPermScores_1e+05Perm_cov_simWithNoise05_Cache", i), {
#         tmp
#     },
#                  cacheDir = paste0(subCache, "/", "rsPermScores_1e+05Perm_cov_simWithNoise05/"), 
#                  recreate=TRUE)
# }

# get histogram of p-values when sampling 300 permutations

nullDistList = COCOA::convertToFromNullDist(a$permRSScores)
gPValDF <- getGammaPVal(rsScores = smallChange[, signalCol], 
                        nullDistList = nullDistList, 
                        signalCol = signalCol, 
                        method = "mme", 
                        realScoreInDist = TRUE,
                        force=FALSE)
hist(nullDistList[[26]]$PC2)
# gPValDF <- apply(X = gPValDF, MARGIN = 2, 
#                 FUN = function(x) p.adjust(p = x, method = correctionMethod))
gPValDF <- cbind(gPValDF, 
                 rsScores[, colnames(rsScores)[!(colnames(rsScores) 
                                                 %in% colsToAnnotate)]])
set.seed(1000)
crossSampleNumber = 500000
gPValList = lapply(X = 1:crossSampleNumber, function(x) x)
pValInd = 17:22
subList = function(myList, innerInd) {
    return(myList[innerInd, ])
}

for (i in 1:crossSampleNumber) {
    randInd = sample(x = 1:100000, size = 300, replace = FALSE)
    myDistList = lapply(X = nullDistList[pValInd], FUN = function(x) subList(x, randInd))
    gPValList[[i]] <- getGammaPVal(rsScores = smallChange[pValInd, signalCol, drop=FALSE], 
                                   nullDistList = myDistList, 
                                   signalCol = signalCol, 
                                   method = "mme", 
                                   realScoreInDist = TRUE,
                                   force=FALSE)
    gPValList[[i]]$rsID = 1:nrow(gPValList[[i]])
}
sampleP = rbindlist(gPValList[1:500000])
# hist(filter(sampleP, rsID == 6)$PC2)
simpleCache(paste0("gammaPSampling300_", i), {
    sampleP
}, assignToVariable = "sampleP")

# 
sP = pivot_longer(data=sampleP, cols=c("PC1", "PC2"), names_to="PC", values_to = "p_value")
sP$rsID = as.factor(sP$rsID)
empPVal = as.data.frame(a$empiricalPVals)[pValInd, signalCol]
empPVal$rsID = as.factor(1:6)
empPVal = pivot_longer(empPVal, cols=c("PC1", "PC2"), names_to = "PC", values_to = "p_value")
tmp = ggplot(data = sP, mapping = aes(x=rsID, y=-log10(p_value))) + geom_boxplot(outlier.shape = NA) + facet_wrap(~PC) +
    geom_crossbar(data=empPVal, 
                  aes(ymin=-log10(p_value), ymax=-log10(p_value)), 
                  fatten=0, color="red") + xlab("Region set") + ylab("Log10 p-value")
tmp

ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "gammaVsEmp300_", .analysisID, ".svg"), 
       plot = tmp, device = "svg", width = plotWidth * .8, height = plotHeight * .8, units = plotUnits)

######
# try sampling 1000 instead of 300
.analysisID = paste0(dataID, "_gauss05")
set.seed(1000)
crossSampleNumber = 500000
gPValList = lapply(X = 1:crossSampleNumber, function(x) x)
subSmallChange = smallChange[pValInd, signalCol, drop=FALSE]
kSamp = replicate(crossSampleNumber, expr = sample(x = 1:100000, 
                                               size = 1000, replace = FALSE))
subDistList = nullDistList[pValInd]
tmpFun = function(randInd) {
    myDistList = lapply(X = subDistList, FUN = function(x) subList(x, randInd))
    getGammaPVal(rsScores = subSmallChange, 
                                   nullDistList = myDistList, 
                                   signalCol = signalCol, 
                                   method = "mme", 
                                   realScoreInDist = TRUE,
                                   force=FALSE)
}

gammaPValsAll = apply(X = kSamp, MARGIN = 2, FUN = tmpFun)

sampleP = rbindlist(gammaPValsAll)
sampleP$rsID = rep(1:6, nrow(sampleP) / length(pValInd))
simpleCache(paste0("gammaPSampling1000_500000"), {
    sampleP
}, assignToVariable = "sampleP")

sP = pivot_longer(data=sampleP, cols=c("PC1", "PC2"), names_to="PC", values_to = "p_value")
sP$rsID = as.factor(sP$rsID)
empPVal = as.data.frame(a$empiricalPVals)[pValInd, signalCol]
empPVal$rsID = as.factor(1:6)
empPVal = pivot_longer(empPVal, cols=c("PC1", "PC2"), names_to = "PC", values_to = "p_value")
tmp = ggplot(data = sP, mapping = aes(x=rsID, y=-log10(p_value))) + geom_boxplot(outlier.shape = NA) + facet_wrap(~PC) +
    geom_crossbar(data=empPVal, 
                  aes(ymin=-log10(p_value), ymax=-log10(p_value)), 
                  fatten=0, color="red") + xlab("Region set") + ylab("Log10 p-value")
tmp
ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "gammaVsEmp1000_", .analysisID, ".svg"), 
       plot = tmp, device = "svg", width = plotWidth * .8, height = plotHeight * .8, units = plotUnits)


##############
# #bumphunter and LOLA
sampleStatus = as.numeric(grepl(pattern = "tumor", 
                                x = colnames(bothHigh), ignore.case = TRUE))
designMat = model.matrix(object = ~ sampleStatus)

simpleCache(paste0("tumorDMR_", dataID, "_gauss025"), {
    tumorDMRL =bumphunter(object = bothLow, design=designMat, chr=COCOA:::grToDt(signalCoord)$chr,
                          pos=COCOA:::grToDt(signalCoord)$start, coef=2, cutoff = 0.05, B=1000)# , type="Beta")
    tumorDMRL
}, assignToVariable = "tumorDMRL")
sum(tumorDMRL$table$fwer <= 0.05)

simpleCache(paste0("tumorDMR_", dataID, "_gauss05"), {
    tumorDMRH =bumphunter(object = bothHigh, design=designMat, chr=COCOA:::grToDt(signalCoord)$chr,
                          pos=COCOA:::grToDt(signalCoord)$start, coef=2, cutoff = 0.05, B=1000)# , type="Beta")
    tumorDMRH
}, assignToVariable = "tumorDMRH")
sum(tumorDMRH$table$fwer <= 0.05)

simpleCache(paste0("lolaResults_", dataID, "_gauss025"), {
    lowResults = dmrLOLA(genomicSignal=bothLow, signalCoord=signalCoord,
                         GRList=GRList, rsAnno=rsAnno, targetVar=sampleStatus,
                         dataID=paste0(dataID, "_gauss025"))
    lowResults
}, assignToVariable = "lResults")

simpleCache(paste0("lolaResults_", dataID, "_gauss05"), {
    hResults = dmrLOLA(genomicSignal=bothHigh, signalCoord=signalCoord,
                         GRList=GRList, rsAnno=rsAnno, targetVar=sampleStatus,
                         dataID=paste0(dataID, "_gauss05"))
    hResults
}, assignToVariable = "hResults")

# COCOA
setLapplyAlias(6)
simpleCache(paste0("cocoaRes_", dataID, "_gauss025_absValT"), {
    smallChange = runCOCOA(genomicSignal = bothLow, signalCoord = signalCoord, GRList = GRList, 
                           signalCol = c("PC1", "PC2"), targetVar = bothLPCA, 
                           variationMetric = "cov", scoringMetric = "regionMean", 
                           absVal = TRUE, centerGenomicSignal = TRUE, centerTargetVar = TRUE)
    smallChange = cbind(smallChange, rsName)
    smallChange
}, assignToVariable = "rsScores025", recreate = FALSE)

simpleCache(paste0("cocoaRes_", dataID, "_gauss05_absValT"), {
smallChange = runCOCOA(genomicSignal = bothHigh, signalCoord = signalCoord, GRList = GRList, 
                       signalCol = c("PC1", "PC2"), targetVar = bothHPCA, 
                       variationMetric = "cov", scoringMetric = "regionMean", 
                       absVal = TRUE, centerGenomicSignal = TRUE, centerTargetVar = TRUE)
smallChange = cbind(smallChange, rsName)
smallChange
}, assignToVariable = "rsScores05", recreate=FALSE)
###################################
# make plots comparing LOLA and COCOA

# identify region sets that overlap with region set of interest
anyOverlap = function(gr1, gr2) {
    fo = findOverlaps(gr1, gr2)
    if (length(subjectHits(fo)) == 0) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

GLL = getOLRegions(GRList = GRList, intGR = signalCoord)
GRList = GLL$GRList


belowCovInd = rsScores025$regionSetCoverage < 100
GRList = GRList[!belowCovInd]
rsScores025 = rsScores025[!belowCovInd, ]
rsScores05 = rsScores05[!belowCovInd, ]
rsName = rsName[!belowCovInd]
rsDescription = rsDescription[!belowCovInd]
lResults = lResults[lResults$filename %in% rsName]

myGR = GRList[["Human_MDA-MB-231-Cells_ESR1,-DBDmut_E2-45min_Katzenellenbogen.bed"]]

# jacOL = sapply(X = GRList, FUN = function(x) getJaccard(rs1 = x, 
#                                                         rs2 = myGR))
# jacOL[jacOL > 1] = 0
# hist(jacOL, breaks=seq(0, 1, 0.01))
# head(names(sort(jacOL, decreasing = TRUE))[2:21])

propOL = sapply(X = GRList, FUN = function(x) percentOL(query = myGR, subject = x))

whichOL = sapply(X = GRList, 
                 FUN = function(x) anyOverlap(x, myGR))
olRSNum = sum(whichOL)
similarRS = names(sort(propOL, decreasing = TRUE))[1:round((olRSNum-1)/100)]
similarRS = similarRS[similarRS != "Human_MDA-MB-231-Cells_ESR1,-DBDmut_E2-45min_Katzenellenbogen.bed"]
olPattern = paste0(similarRS, collapse = "|")


loopNames = c("025", "05")

for (j in seq_along(loopNames)) {
    
    rsScores = get(paste0("rsScores", loopNames[j]))
    .analysisID = paste0(dataID, "_gauss", loopNames[j], "_absValT")
    
    # region score distribution
    for (i in paste0("PC", 1:2)) {
        
        a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = i, 
                              pattern = c("Human_MDA-MB-231-Cells_ESR1,-DBDmut_E2-45min_Katzenellenbogen.bed", 
                                          olPattern), 
                              patternName = c("Region set of interest", "Some overlap with region set")) +
            theme(legend.position = c(0.15, 0.15)) +
            scale_color_manual(values = c("gray", "blue", "red")) + 
            xlab(paste0("Region set rank (", i, ")")) + 
            theme(axis.title.y = element_blank(), 
                  legend.text = element_blank(), 
                  legend.title = element_blank(), 
                  legend.position = "none") +
            scale_x_continuous(breaks = c(0, 1000, 2000), 
                               labels= c("0", "1000", "2000"), limits=c(-25, nrow(rsScores) + 25))
        
        a 
        ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", i, "_", .analysisID, ".svg"), 
               plot = a, device = "svg", width = plotWidth * .8, height = plotHeight * .8, units = plotUnits)
    }
}


# making a plot with legend
i="PC1"
a = plotAnnoScoreDist(rsScores = rsScores, colsToPlot = i, 
                      pattern = c("Human_MDA-MB-231-Cells_ESR1,-DBDmut_E2-45min_Katzenellenbogen.bed", olPattern), 
                      patternName = c("Region set of interest", "Some overlap with region set")) +
    theme(legend.position = c(0.15, 0.15)) +
    scale_color_manual(values = c("gray", "blue", "red")) + 
    xlab(paste0("Region set rank (", i, ")")) 

a 
ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", i, "_", .analysisID, "_withLegend.svg"), 
       plot = a, device = "svg", width = plotWidth, height = plotHeight, units = plotUnits)

#################
# plots for LOLA results
.analysisID = paste0(dataID, "_gauss025")

lResults = arrange(lResults, desc(oddsRatio))
lResults$oddsRank = 1:nrow(lResults)
rsGroup = rep("Other", nrow(lResults))
rsGroup[lResults$filename == "Human_MDA-MB-231-Cells_ESR1,-DBDmut_E2-45min_Katzenellenbogen.bed"] = "Region set of interest"
rsGroup[lResults$filename %in% similarRS] = "Some overlap with region set"
lResults$rsGroup = factor(rsGroup)

myPlot = ggplot(data = lResults, aes(x=oddsRank, y = oddsRatio)) + geom_point(shape=3, aes(col=rsGroup), alpha=0.5) +
    scale_color_manual(values = c("gray", "blue", "red")) + xlab(paste0("Region set rank")) + 
    theme(legend.text = element_blank(), 
          legend.title = element_blank(), 
          legend.position = "none") +
    scale_x_continuous(breaks = c(0, 1000, 2000), 
                       labels= c("0", "1000", "2000"), limits=c(-25, nrow(rsScores) + 25)) +
    ylab("Odds ratio")

myPlot

ggsave(filename = paste0(Sys.getenv("PLOTS"), plotSubdir, "annoScoreDist_", "LOLA", "_", .analysisID, ".svg"), 
       plot = myPlot, device = "svg", width = plotWidth * .8, height = plotHeight * .8, units = plotUnits)
#######################################################################


##########################################################################
