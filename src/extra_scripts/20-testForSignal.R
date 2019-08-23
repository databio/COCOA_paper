# see whether there is any variation in regions of interest
# if there is no signal there in the first place, we could not use those
# region sets as true positives for testing the methods
library(COCOA)

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
plotSubDir = "20testSignal/"

if (!dir.exists(ffPlot(plotSubDir))) {
    dir.create(ffPlot(plotSubDir))
}



###############################################################################

simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList")

###########################################################
# reading in the region sets
# load LOLA database
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))
names(GRList) = rsName

#################################################################

dataID = "657" # 657 patients with both ER and PGR info in metadata, 692 total
allMPCAString = "allMPCA_657" #  "allMPCA_657"
top10MPCAString = "top10MPCA_657"

simpleCache("rsEnrichmentRSMean_657", assignToVariable = "rsScores")

# get top region sets for DNA methylation results
topScoreResultsInd = rsRankingIndex(rsScores, "PC1")$PC1[1:20]
topGRList = GRList[rsScores$rsName[topScoreResultsInd]]

# loading PCA and combining components that could separate ER=/-
# for rsEnrichment, PCs 1 and 4 could separate ER+/-
simpleCache(allMPCAString, assignToVariable = "mPCA")
brcaLoadings = mPCA$rotation
brcaCoord = brcaMList[["coordinates"]]

############################################################################
# protein data

# proteins of interest (we have a region set related to them):
# NFKB, GATA3, GATA6, ERALPHA_pS118, AR, CJUN_pS73, CMYC, , 
# BETACATENIN, PR, SMAD1, SMAD3, STAT3_pY705, STAT5ALPHA, ETS1, FOXM1, IRF1, 
# interesting but not necessarily clear what the "right" answer is: 
# CD31, BRAF, CYCLIND1, RICTOR, CD49B, CD20
mainP = c("NFKBP65_pS536", "GATA3", "GATA6", "ERALPHA_pS118", "ERALPHA", "AR", 
          "CJUN_pS73", "CMYC", "BETACATENIN", "PR", "SMAD1", "SMAD3", "STAT3_pY705", 
          "STAT5ALPHA", "ETS1", "FOXM1", "IRF1")
alsoOfInterest = c("CD31", "BRAF", "CYCLIND1", "RICTOR", "CD49B", "CD20")
proteinOfInterest = c(mainP, alsoOfInterest)

prot = read.csv(file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/TCGA-BRCA-L4.csv"))
prot$subject_ID = substr(prot$Sample_ID, start = 1, stop = 12)

# figuring out which patients are represented more than once
idt = table(prot$subject_ID)
# head(sort(idt, decreasing = TRUE))
multiInd = which(idt > 1)
multiPatient = names(multiInd)

# remove patients with more than one sample
prot = prot[!(prot$subject_ID %in% multiPatient), ]
row.names(prot) = prot$subject_ID

sum(colnames(brcaMList[[2]]) %in% prot$subject_ID)
sharedSamples = colnames(brcaMList[[2]])[colnames(brcaMList[[2]]) %in% prot$subject_ID]

filteredMData = brcaMList[[2]][, sharedSamples]
filtProt = prot[sharedSamples, ]

#### convert DNA methylation matrix to correlation matrix
# calculate correlation
featurePCCor = createCorFeatureMat(dataMat = filteredMData, 
                                   featureMat = as.matrix(filtProt[, proteinOfInterest]), 
                                   centerDataMat=TRUE, centerFeatureMat=TRUE)
colnames(featurePCCor) <- proteinOfInterest
hist(featurePCCor[, "GATA3"])
############################################################################
# does DNAm in region sets for protein show correlation with protein expression

if (!dir.exists(ffPlot(paste0(plotSubDir, "QQPlots/")))) {
    dir.create(ffPlot(paste0(plotSubDir, "QQPlots/")))
}

# get background distribution for all CpGs
grDevices::pdf(file = ffPlot(paste0(plotSubDir, "allCpG_DNAmCorHist.pdf")))
for (i in seq_along(proteinOfInterest)) {
    hist(abs(featurePCCor[, proteinOfInterest[i]]), 
         breaks = seq(0, 1, 0.025), main = proteinOfInterest[i])
}
dev.off()

pValList = list()

colMeans(abs(featurePCCor))
corByRegionList = list()

proteinOfInterest = c("NFKBP65_pS536", "GATA3", "GATA6", "ERALPHA_pS118", "ERALPHA", 
                      "CJUN_pS73", "CMYC", "BETACATENIN", "PR", "SMAD1", "SMAD3", "STAT3_pY705", 
                      "STAT5ALPHA", "ETS1", "FOXM1", "IRF1")

protSearchTerm = c("NFKB|nfkappa", "GATA3", "GATA6", "ERALPHA|eraa|esr1", "ERALPHA|eraa|esr1", 
                   "JUN", "MYC", "tcf3|tcfl7", "pgr", "SMAD1", "SMAD3", "STAT3", 
                   "STAT5", "ETS1", "FOXM1", "IRF1")

for (i in seq_along(proteinOfInterest)) {
    thisGRList = GRList[grep(pattern = protSearchTerm[i], x = names(GRList), ignore.case = TRUE)]
    PCsToAnnotate = proteinOfInterest[i]
    loadingMat = as.matrix(featurePCCor[, proteinOfInterest[i]])
    colnames(loadingMat) = proteinOfInterest[i]
    corByRegionList[[i]] = lapply(X = thisGRList, FUN = function(x) COCOA:::averageByRegion(loadingMat = loadingMat, 
                                                                     signalCoord = brcaCoord, 
                                                                     regionSet = x, PCsToAnnotate = PCsToAnnotate))
    
   
        
    grDevices::pdf(file = ffPlot(paste0(plotSubDir, proteinOfInterest[i], "_DNAmCorHist.pdf")))
    for (j in seq_along(corByRegionList[[i]])) {
        hist(as.numeric(as.matrix(corByRegionList[[i]][[j]][, proteinOfInterest[i], with=FALSE])), 
                   breaks = seq(0, 1, 0.025), main = proteinOfInterest[i])
    }
    dev.off()
    
    # quantile quantile plots
    grDevices::pdf(file = ffPlot(paste0(plotSubDir, "QQPlots/", proteinOfInterest[i], "_DNAmCorQQPlot.pdf")))
    for (j in seq_along(corByRegionList[[i]])) {
        qqplot(x = abs(featurePCCor[, proteinOfInterest[i]]), 
               y = as.numeric(as.matrix(corByRegionList[[i]][[j]][, proteinOfInterest[i], with=FALSE])), 
               xlab = "reference", ylab = proteinOfInterest[i], main= "QQ plot of correlation")
        lines(x = seq(0, 1, 0.5), y = seq(0, 1, 0.5))
    }
    dev.off()
    pValList[[i]] = list()
    for (j in seq_along(corByRegionList[[i]])) {
        pValList[[i]][[j]] = wilcox.test(x = abs(featurePCCor[, proteinOfInterest[i]]), 
                                         y= as.numeric(as.matrix(corByRegionList[[i]][[j]][, proteinOfInterest[i], with=FALSE])), 
                                         conf.int = TRUE)
    }

    
        
    # grDevices::pdf(file = ffPlot(paste0(plotSubDir, proteinOfInterest[i], "_DNAmCor.pdf")), 
    #                width = 11, 
    #                height = 8.5 * length(corByRegionList[[i]]))
    # grid.newpage()
    # for (j in seq_along(corByRegionList[[i]])) {
    #     multiHM <- grid.grabExpr(draw(Heatmap(matrix = as.matrix(corByRegionList[[i]][[j]][, proteinOfInterest[i], with=FALSE]), 
    #                                                 col = c("white", "black"))))
    #     
    #     pushViewport(viewport(y = unit((8.5 * length(corByRegionList[[i]])) -
    #                                        (j - 1) * 8.5, "in"), 
    #                           height = unit(8, "in"), just = "top"))
    #     grid.draw(multiHM)
    #     popViewport()
    # }
    # dev.off()
    
}



#############################################################################
# testing simple variance

cpgVar = apply(t(brcaMList$methylProp), MARGIN = 2, FUN = var)
length(cpgVar)
PCsToAnnotate = "cpgVar"

stat3RSInd = grep(pattern = "stat3", x = names(GRList), ignore.case = TRUE)
nfkbRSInd = grep(pattern = "nfkb|nfkappa", x = names(GRList), ignore.case = TRUE)

varByRegion = COCOA:::averageByRegion(loadingMat = data.frame(cpgVar = cpgVar), signalCoord = brcaCoord, regionSet = GRList[[stat3RSInd[2]]], PCsToAnnotate = PCsToAnnotate)

Heatmap(matrix = as.matrix(varByRegion$cpgVar), col = c("white", "red"))

# does PCA of just regions for region set of interest return good results
rsPCA = dimRedOnRS(regionSet = GRList[[stat3RSInd[2]]], methylData = filteredMData, mCoord = brcaCoord, drMethod = "pca")

protOfInt = "NFKBP65_pS536"
protOfInt = "STAT3_pY705"
colorClusterPlots(clusteredDF = cbind(rsPCA$x, STAT3_pY705=filtProt[, protOfInt]), 
                  plotCols = c("PC1", "PC2"), colorByCols = protOfInt, alphaVal = 1)
