# library(projectInit)

# project.init(codeRoot = paste0(Sys.getenv("CODE"), "COCOA/R/"), dataDir = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/"))
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
# library(fastICA)

# 
setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)
plotSubdir = "02_COCOA_DNAm/"

# DNA methylation data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))

#############################################################################



simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList")

# simpleCache("mPCA_692", {
#     mPCA_all_samples = prcomp(t(brcaMList$methylProp))
# })


#####
# 
# foxa1 = GRList1[[grep(pattern = "FoxA1_E2-45min_Brown", 
#                       x = lolaCoreRegionAnno$filename)]]
# topFoxa1 = getTopRegions(loadingMat = allMPCA_657$rotation,
#                          signalCoord = COCOA:::dtToGr(brcaMList$coordinates), 
#                          regionSet = foxa1, 
#                          PCsToAnnotate = c("PC1", "PC4"), returnQuantile = TRUE)
# gata3 = GRList1[[grep(lolaCoreRegionAnno$filename, pattern = "SydhMcf7Gata3Ucd")]]
# topGata3 = getTopRegions(loadingMat = allMPCA_657$rotation,
#                          signalCoord = COCOA:::dtToGr(brcaMList$coordinates), 
#                          regionSet = gata3, 
#                          PCsToAnnotate = c("PC1", "PC4"), returnQuantile = TRUE)
# er = GRList1[[grep(lolaCoreRegionAnno$filename, pattern = "ESR1_E2-45min_Brown")]]
# topER = getTopRegions(loadingMat = allMPCA_657$rotation,
#                          signalCoord = COCOA:::dtToGr(brcaMList$coordinates), 
#                          regionSet = er, 
#                          PCsToAnnotate = c("PC1", "PC4"), returnQuantile = TRUE)
# h3r17 = GRList1[[grep(lolaCoreRegionAnno$filename, pattern = "H3R17me2_E2-45min_Brown")]]
# topH3R17 = getTopRegions(loadingMat = allMPCA_657$rotation,
#                       signalCoord = COCOA:::dtToGr(brcaMList$coordinates), 
#                       regionSet = h3r17, 
#                       PCsToAnnotate = c("PC1", "PC4"), returnQuantile = TRUE)
# ezh2 = GRList1[[grep(lolaCoreRegionAnno$filename, pattern = "HmecEzh239875")]]
# topEZH2 = getTopRegions(loadingMat = allMPCA_657$rotation,
#                          signalCoord = COCOA:::dtToGr(brcaMList$coordinates), 
#                          regionSet = ezh2, 
#                          PCsToAnnotate = c("PC1", "PC4"), returnQuantile = TRUE)
# topRS = c(topFoxa1, topGata3, topER, topH3R17, topEZH2)
# topRSCombined = reduce(topRS)
# save(topRSCombined, file="topRSCombined.RData")
########

#restrict patients included in this analysis
patientMetadata = patientMetadata[patientMetadata$subject_ID %in% 
                                      colnames(brcaMList[["methylProp"]]), ]


# patientMetadata should have already screened out patients without ER/PGR status
# resulting in 657 patients
hasER_PGR_IDs = patientMetadata[, subject_ID]
filteredMData = brcaMList[["methylProp"]][, 
                                          colnames(brcaMList[["methylProp"]]) %in% hasER_PGR_IDs] 

###########################################################
# reading in the region sets
# load LOLA database
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))

#################################################################

dataID = "657" # 657 patients with both ER and PGR info in metadata, 692 total
allMPCAString = "allMPCA_657" #  "allMPCA_657"
top10MPCAString = "top10MPCA_657"
rsName = c("GSM2305313_MCF7_E2_peaks_hg38.bed", 
           lolaCoreRegionAnno$filename,
           roadmapRegionAnno$filename,
           motifRegionAnno$filename)
rsDescription = c("ER ChIPseq, MCF7 cell line, estradiol stimulation",
                  lolaCoreRegionAnno$description,
                  roadmapRegionAnno$description,
                  motifRegionAnno$description)
           
# loading PCA and combining components that could separate ER=/-
# for rsEnrichment, PCs 1 and 4 could separate ER+/-
simpleCache(allMPCAString, assignToVariable = allMPCAString)
simpleCache(top10MPCAString, assignToVariable = "top10MPCA")
# getting new axis in direction that separates, normalizing to length one
PC1m4 = (1/sqrt(2)) * allMPCA$rotation[, "PC1"] - (1/sqrt(2)) * allMPCA$rotation[, "PC4"]
allMPCA$rotation = cbind(allMPCA$rotation, PC1m4)
# also adding combination that could separate ER+/- in top10MPCA which should not here
PC1p3 = (1/sqrt(2)) * allMPCA$rotation[, "PC1"] + (1/sqrt(2)) * allMPCA$rotation[, "PC3"]
allMPCA$rotation = cbind(allMPCA$rotation, PC1p3)

# for rsEnrichmentTop10, PCs 1 and 3 could separate ER+/-
PC1p3 = (1/sqrt(2)) * top10MPCA$rotation[, "PC1"] + (1/sqrt(2)) * top10MPCA$rotation[, "PC3"]
top10MPCA$rotation = cbind(top10MPCA$rotation, PC1p3)
# also adding combination that could separate ER+/- in allMPCA which should not here
PC1m4 = (1/sqrt(2)) * top10MPCA$rotation[, "PC1"] - (1/sqrt(2)) * top10MPCA$rotation[, "PC4"]
top10MPCA$rotation = cbind(top10MPCA$rotation, PC1m4)
# make COCOA_pipeline be able to take PCA object? otherwise create new cache with PC values
# specialPCEnr = COCOA_pipeline(mData=filteredMData, coordinates=brcaMList[["coordinates"]], 
#                               GRList=GRList, useCache=TRUE, 
#                               allMPCAString=allMPCAString, top10MPCAString = top10MPCAString, 
#                               rsName = rsName, rsDescription = rsDescription)



# gives output of rsEnrichment from PCA of all shared cytosines
# and rsEnrichmentTop10 from PCA of 10% most variable shared cytosines
source(paste0(Sys.getenv("CODE"),"/aml_e3999/src/COCOA_pipeline.R"))
# brcaMList = filteredMData
enrichResults = COCOA_pipeline(mData=NULL, coordinates=brcaMList[["coordinates"]], 
               GRList=GRList, 
               PCsToAnnotate =paste0("PC", 1:10), #  c("PC1m4", "PC1p3", 
               scoringMetric = "rsMean",
               pcaCache=FALSE, 
               allMPCACacheName=allMPCAString,  
               overwritePCACaches = FALSE, 
               allMPCA = get(allMPCAString), 
               rsName = rsName, rsDescription = rsDescription,
               rsEnCache = TRUE, rsEnCacheName = ,
               
               overwriteResultsCaches = TRUE) 

rsEnrichment = enrichResults[[1]]
write.csv(x = rsEnrichment, 
              file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/PC_Enrichment_All_Shared_Cs_692_rsMean.csv"),
              quote = FALSE, row.names = FALSE)

# running again for top 10 most variable CpGs
enrichResults = COCOA_pipeline(mData=NULL, coordinates=NULL, 
                               GRList=GRList, 
                               PCsToAnnotate =paste0("PC", 1:10), #  c("PC1m4", "PC1p3", 
                               scoringMetric = "rsMean",
                               pcaCache=FALSE, 
                               allMPCACacheName=top10MPCAString, 
                               overwritePCACaches = FALSE, 
                               allMPCA = top10MPCA,
                               rsName = rsName, rsDescription = rsDescription,
                               rsEnCache = TRUE, rsEnCacheName = "rsEnrichmentTop10_657",
                               overwriteResultsCaches = TRUE) 

rsEnrichmentTop10 = enrichResults[[2]]
write.csv(x = rsEnrichmentTop10, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/PC_Enrichment_Top_10%_Variable_Cs_657.csv"),
          quote = FALSE, row.names = FALSE)


#################################################################
# run COCOA again but with wilcoxon rank sum scoring method

# gives output of rsEnrichment from PCA of all shared cytosines
# and rsEnrichmentTop10 from PCA of 10% most variable shared cytosines
enrichResults = COCOA_pipeline(mData=filteredMData, coordinates=brcaMList[["coordinates"]], 
                               GRList=GRList, 
                               PCsToAnnotate = c(paste0("PC", 1:4)), 
                               scoringMetric = "rankSum",
                               pcaCache=FALSE, 
                               allMPCACacheName=allMPCAString,  
                               overwritePCACaches = FALSE, 
                               allMPCA = allMPCA, 
                               rsName = rsName, rsDescription = rsDescription,
                               rsEnCache = TRUE, rsEnCacheName = "rsEnrichment_ranksum_657",
                               overwriteResultsCaches = TRUE) 

rsEnrichment = enrichResults
write.csv(x = rsEnrichment, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/PC_Enrichment_All_Shared_Cs_ranksum_657.csv"),
          quote = FALSE, row.names = FALSE)

# for top 10% variable CpGs
enrichResults = COCOA_pipeline(mData=filteredMData, coordinates=brcaMList[["coordinates"]], 
                               GRList=GRList, 
                               PCsToAnnotate = c(paste0("PC", 1:4)), 
                               scoringMetric = "rankSum",
                               pcaCache=FALSE, 
                               allMPCACacheName=top10MPCAString,  
                               overwritePCACaches = FALSE, 
                               allMPCA = top10MPCA,
                               rsName = rsName, rsDescription = rsDescription,
                               rsEnCache = TRUE, rsEnCacheName = "rsEnrichmentTop10_ranksum_657",
                               overwriteResultsCaches = TRUE) 
rsEnrichmentTop10 = enrichResults
write.csv(x = rsEnrichmentTop10, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/PC_Enrichment_Top_10%_Variable_Cs_ranksum_657.csv"),
          quote = FALSE, row.names = FALSE)

#################################################################
# run COCOA again but with mean of individual CpGs instead of mean of regions

# gives output of rsEnrichment from PCA of all shared cytosines
# and rsEnrichmentTop10 from PCA of 10% most variable shared cytosines
enrichResults = COCOA_pipeline(mData=filteredMData, coordinates=brcaMList[["coordinates"]], 
                               GRList=GRList, 
                               PCsToAnnotate = c(paste0("PC", 1:6)), 
                               scoringMetric = "cpgMean",
                               pcaCache=FALSE, 
                               allMPCACacheName=allMPCAString,  
                               overwritePCACaches = FALSE, 
                               allMPCA = allMPCA, 
                               rsName = rsName, rsDescription = rsDescription,
                               rsEnCache = TRUE, rsEnCacheName = "rsEnrichment_cpgMean",
                               overwriteResultsCaches = TRUE) 

rsEnrichment = enrichResults



write.csv(x = rsEnrichment, 
          file = paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/sheets/PC_Enrichment_All_Shared_Cs_rawCpG_657.csv"),
          quote = FALSE, row.names = FALSE)

enrichResults = COCOA_pipeline(mData=filteredMData, coordinates=brcaMList[["coordinates"]], 
                               GRList=GRList, 
                               PCsToAnnotate = c(paste0("PC", 1:6)), 
                               scoringMetric = "cpgMean",
                               pcaCache=FALSE, 
                               allMPCACacheName=top10MPCAString,  
                               overwritePCACaches = FALSE, 
                               allMPCA = top10MPCA, 
                               rsName = rsName, rsDescription = rsDescription,
                               rsEnCache = TRUE, rsEnCacheName = "rsEnrichmentTop10_cpgMean",
                               overwriteResultsCaches = TRUE) 
rsEnrichmentTop10 = enrichResults
write.csv(x = rsEnrichmentTop10, 
          file = dirData("analysis/sheets/PC_Enrichment_Top_10%_Variable_Cs_cpgMean_657.csv"),
          quote = FALSE, row.names = FALSE)

######################################################
# visualizing PCA

# load PCA data
simpleCache(allMPCAString, assignToVariable = "allMPCA")
simpleCache(top10MPCAString, assignToVariable = "top10MPCA")

# shorten metadata for visualization
patientMetadata$menopause_status = sub(pattern = "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)", 
            "Pre", x = patientMetadata$menopause_status, fixed = TRUE)
patientMetadata$menopause_status = sub(pattern = "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)", 
                                       "Post", x = patientMetadata$menopause_status, fixed = TRUE)
patientMetadata$menopause_status = sub(pattern = "Indeterminate (neither Pre or Postmenopausal)", 
                                       "Indeterminate", x = patientMetadata$menopause_status, fixed = TRUE)
patientMetadata$menopause_status = sub(pattern = "Peri (6-12 months since last menstrual period)", 
                                       "Peri", x = patientMetadata$menopause_status, fixed = TRUE)
patientMetadata$race = sub(pattern = "BLACK OR AFRICAN AMERICAN", 
                                       "BLACK", x = patientMetadata$race, fixed = TRUE)
patientMetadata$race = sub(pattern = "AMERICAN INDIAN OR ALASKA NATIVE", 
                           "NATIVE AM. OR ALASKAN", x = patientMetadata$race, fixed = TRUE)

# getting PC scores manually so artificial PCs will be included (PC1m4 and PC1p3)
centeredPCAMeth = t(apply(filteredMData, 2, function(x) x - allMPCA$center)) # center first 
reducedValsPCA = centeredPCAMeth %*% allMPCA$rotation
pcaValDF = as.data.frame(reducedValsPCA)

# add annotation information
pcaWithAnno = cbind(pcaValDF, patientMetadata)

colorByCols = colnames(patientMetadata)[!(colnames(patientMetadata) %in% "subject_ID")]

if (!dir.exists(paste0(Sys.getenv("PLOTS"), "/allMPCA_PCA_Plots"))) {
    dir.create(paste0(Sys.getenv("PLOTS"), "/allMPCA_PCA_Plots"), recursive = TRUE)
}

PCsToPlot = c("PC1m4", "PC1p3", paste0("PC", 1:10))
for (i in seq_along(PCsToPlot)) {
    for (j in seq_along(PCsToPlot)) {
        multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                               plotCols = c(PCsToPlot[i], PCsToPlot[j]), 
                                               colorByCols=colorByCols)
        ggplot2::ggsave(filename=paste0(Sys.getenv("PLOTS"), "/allMPCA_PCA_Plots/multiColorPCAPlots_allMPCA_", 
                                        PCsToPlot[i], "x", PCsToPlot[j], 
                                        ".pdf"), plot = multiColorPCAPlots, device = "pdf",
                        limitsize=FALSE)
        
    }
    
}


# for PCA of top 10% variable cytosines
# add annotation information
pcaWithAnno = cbind(as.data.table(top10MPCA$x), patientMetadata)

colorByCols = colnames(patientMetadata)[!(colnames(patientMetadata) %in% "subject_ID")]
for (i in 2:6) {
    multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                           plotCols = c("PC1", paste0("PC", i)), 
                                           colorByCols=colorByCols)
    ggplot2::ggsave(filename=paste0(Sys.getenv("PLOTS"), paste0("multiColorPCAPlots_", top10MPCAString, "_1", i), 
                                    ".pdf"), plot = multiColorPCAPlots, device = "pdf",
                    limitsize=FALSE)
}
for (i in c(2, 4:6)) {
    multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                           plotCols = c("PC3", paste0("PC", i)), 
                                           colorByCols=colorByCols)
    ggplot2::ggsave(filename=paste0(Sys.getenv("PLOTS"), paste0("multiColorPCAPlots_", top10MPCAString, "_3", i), 
                                    ".pdf"), plot = multiColorPCAPlots, device = "pdf",
                    limitsize=FALSE)
}


###################################################################################
# # visualizing ICA
# # treating cytosines as observations and patients as dimensions 
# allMICA = fastICA(X = brcaMList$methylProp, n.comp = 5)#, alg.typ = "deflation")
# 
# i=2
# cpgToPlotNum=20000
# pdf(paste0(Sys.getenv("PLOTS"), "ICA_cpg_plots.pdf"))
# for (i in 1:ncol(allMICA$S)) {
#     for (j in 1:ncol(allMICA$S)) {
#         cpgToPlot = sample(1:nrow(allMICA$S), cpgToPlotNum)
#         plot(allMICA$S[cpgToPlot, i], allMICA$S[cpgToPlot, j])
#     }
# }
# dev.off()
# plot(allMICA$S[1:cpgToPlot, i], allMICA$S[1:cpgToPlot, 2])
# plot(allMICA$S[1:cpgToPlot, i], allMICA$S[1:cpgToPlot, 3])
# plot(allMICA$S[1:cpgToPlot, i], allMICA$S[1:cpgToPlot, 4])
# plot(allMICA$S[1:cpgToPlot, i], allMICA$S[1:cpgToPlot, 5])
# 
# # # add annotation information
# # icaWithAnno = cbind(as.data.table(allMICA$x), patientMetadata[dataSplit, ])
# # 
# # colorByCols = colnames(patientMetadata)[!(colnames(patientMetadata) %in% "subject_ID")]
# # for (i in 2:6) {
# #     multiColorICAPlots = colorClusterPlots(icaWithAnno, 
# #                                            plotCols = c("PC1", paste0("PC", i)), 
# #                                            colorByCols=colorByCols)
# #     ggplot2::ggsave(filename=paste0(Sys.getenv("PLOTS"), paste0("multiColorICAPlots1", i), 
# #                                     ".pdf"), plot = multiColorICAPlots, device = "pdf",
# #                     limitsize=FALSE)
# # }
# # for (i in c(2, 4:6)) {
# #     multiColorICAPlots = colorClusterPlots(icaWithAnno, 
# #                                            plotCols = c("PC3", paste0("PC", i)), 
# #                                            colorByCols=colorByCols)
# #     ggplot2::ggsave(filename=paste0(Sys.getenv("PLOTS"), paste0("multiColorICAPlots3", i), 
# #                                     ".pdf"), plot = multiColorICAPlots, device = "pdf",
# #                     limitsize=FALSE)
# # }
# 




# visualizing 


######################################################################################
# analysis of whether certain subsets of region sets are the variable ones

# erSet = GRangesList(MIRA:::dtToGr(erSet))
# genomeLoadings = cbind(coordinates, as.data.table(allMPCA$rotation))
# regionAv = averageByRegion(BSDT = genomeLoadings[, .(chr, start, PC1, PC2, PC3, PC4)], regionsGRL = erSet, 
#                 jCommand = MIRA:::buildJ(c("PC1", "PC2", "PC3", "PC4"), "mean"),
#                 hasCoverage = FALSE)
simpleCache("regionAv101", {
    regionAv = lapply(GRList[1:100], function(x) averageByRegion(loadingMat = allMPCA$rotation, coordinateDT= coordinates, GRList = x, 
                                                                  PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")))
    names(regionAv) <- names(GRList)[1:100]
    regionAv
})
grep(pattern = "A549P300", names(regionAv))
# two region sets that show small peaks in the middle of dips
# these have some regions that have very high loadings which probably drive that
hist(regionAv$wgEncodeAwgTfbsHaibA549P300V0422111Etoh02UniPk.narrowPeak$PC5)
hist(regionAv$wgEncodeAwgTfbsHaibA549Tcf12V0422111Etoh02UniPk.narrowPeak$PC2)
# a region set with a strong peak, has a long tail but very robust
hist(regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`$PC6)
# are the high regions in one PC the same as the high regions in another PC?
View(regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[order(PC6, decreasing = TRUE), ])
plot(regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[, .(PC1, PC3)])
plot(regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[, .(PC1, PC4)])
plot(regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[, .(PC3, PC4)])
cor(x = regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[ PC1 > .002, PC1],
    y = regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[ PC1 > .002, PC3])

# LOLA analysis of regions that are high in one PC compared to another for 
# Human_MCF−7_ESR1_E2−45min_Brown.bed which had good peaks for PC1, 3, 4, 5
plot(regionAv$`Human_MCF-7_ESR1_E2-45min_Brown.bed`[, .(PC1, PC3)])
esr1 = regionAv$`Human_MCF-7_ESR1_E2-45min_Brown.bed`
pc1Reg = esr1[PC1 > 0.001, .(chr, start, end)]
pc2Reg = esr1[PC2 > 0.001, .(chr, start, end)]
pc3Reg = esr1[PC3 > 0.001, .(chr, start, end)]
RGenomeUtils::writeBed(esr1[, .(chr, start, end)], filename = "esr1_all.bed")
RGenomeUtils::writeBed(pc1Reg, filename = "esr1_PC1.bed")
RGenomeUtils::writeBed(pc2Reg, filename = "esr1_PC2.bed")
RGenomeUtils::writeBed(pc3Reg, filename = "esr1_PC3.bed")
fos = regionAv$`Human_MCF-7_c-Fos_E2-45min-3hr_Liu.bed`
pc1Reg = fos[PC1 > 0.001, .(chr, start, end)]
pc2Reg = fos[PC2 > 0.001, .(chr, start, end)]
pc3Reg = fos[PC3 > 0.001, .(chr, start, end)]
RGenomeUtils::writeBed(fos[, .(chr, start, end)], filename = "fos_all.bed")
RGenomeUtils::writeBed(pc1Reg, filename = "fos_PC1.bed")
RGenomeUtils::writeBed(pc2Reg, filename = "fos_PC2.bed")
RGenomeUtils::writeBed(pc3Reg, filename = "fos_PC3.bed")
# findOverlaps(dtToGr(pc1RegEsr), dtToGr(pc1RegFos))

#############################################################################
# generate plots of methylation distribution across patients
# one plot per cpg
set.seed(1234)
numCpgsToPlot = 10000
cpgsToPlot = sample(x = 1:nrow(brcaMList$methylProp), numCpgsToPlot)
histList = list()
# grDevices::pdf(paste0(Sys.getenv("PLOTS"), "cpg_methyl_distr_across_samples.pdf"))
for (i in 1:numCpgsToPlot) {
    histList[[i]] = t(hist(brcaMList$methylProp[cpgsToPlot[i], ], breaks = seq(from=0, to=1, by = 0.1))$count)
}
countMat = do.call(rbind, histList)
# dev.off()
# scale max to 1 for heatmap
histList2 = lapply(histList, FUN = function(x) x / max(x))

countMat2 = do.call(rbind, histList2)

grDevices::pdf(paste0(Sys.getenv("PLOTS"), "cpg_methyl_distr_across_samples_heatmap.pdf"))
Heatmap(countMat, cluster_columns = FALSE, cluster_rows = TRUE)
dev.off()
grDevices::pdf(paste0(Sys.getenv("PLOTS"), "cpg_methyl_distr_across_samples_heatmap_scaled.pdf"))
Heatmap(countMat2, cluster_columns = FALSE, cluster_rows = TRUE)
dev.off()



###################################################################################

# # could also do the analysis for hormone receptors in breast cancer:
# # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4754852/
# load("allBRCAexpression.RData")
# # figuring out which patients are ER+ (inexact estimate based on gene expression data)
# ER1Ind = grep(pattern = "ENSG00000091831", x = myExprDT[, "Gene"],
#      ignore.case = TRUE) # ESR1
# hist(as.numeric(myExprDT[ER1Ind, 2:ncol(myExprDT)]),
#      breaks = seq(0,300,1))
# sum(as.numeric(myExprDT[ER1Ind, 2:ncol(myExprDT)]) < 1) / (ncol(myExprDT)-1)
# # progesterone receptor positive
# PGRInd = grep(pattern = "ENSG00000082175", x = myExprDT[, "Gene"],
#               ignore.case = TRUE) # PGR
# hist(as.numeric(myExprDT[PGRInd, 2:ncol(myExprDT)]),
#      breaks = seq(0,400,1))
# sum(as.numeric(myExprDT[PGRInd, 2:ncol(myExprDT)]) < 2) / (ncol(myExprDT)-1)


# References
# https://www.ncbi.nlm.nih.gov/pubmed/17616709/: ER transcriptional network


## 
# for (i in 1:30) {
#     
# 
# t = apply(X = allMPCA$rotation[, 1:10], 2, FUN = sd)
# }


# testing whether loadings are really correlation of a given feature with that PC score for a given PC
cpgMeans = rowMeans(filteredMData)
# centering before calculating correlation
centeredMData = apply(X = filteredMData, MARGIN = 2, function(x) x - cpgMeans)
pcScores = allMPCA_657$x
pcLoadings = allMPCA_657$rotation
all(colnames(filteredMData) == row.names(pcScores))
featureCorPC1 = apply(X = centeredMData, 1, FUN = function(x) cor(x = x, pcScores[, "PC1"]))
loadCorRatioPC1 = pcLoadings[, "PC1"] / featureCorPC1^2
hist(featureCorPC1)
hist(loadCorRatioPC1[(loadCorRatioPC1 < 1) & (loadCorRatioPC1 > -1)])

featureCorPC2 = apply(X = centeredMData, 1, FUN = function(x) cor(x = x, pcScores[, "PC2"]))
loadCorRatioPC2 = pcLoadings[, "PC2"] / featureCorPC2
hist(loadCorRatioPC2)


################################################################################
# COCOA with PCA done with variance normalization (scale=TRUE)

# patientMetadata should have already screened out patients without ER/PGR status
# resulting in 657 patients
hasER_PGR_IDs = patientMetadata[, subject_ID]
filteredMData = brcaMList[["methylProp"]][, 
                                          colnames(brcaMList[["methylProp"]]) %in% hasER_PGR_IDs] 
patientMetadata = patientMetadata[hasER_PGR_IDs %in% colnames(brcaMList[["methylProp"]])]

normPCA = prcomp(x = t(filteredMData), center = TRUE, scale. = TRUE)

pcaWithAnno = cbind(normPCA$x, patientMetadata)

colorByCols = colnames(patientMetadata)[!(colnames(patientMetadata) %in% "subject_ID")]

if (!dir.exists(ffPlot(paste0(plotSubdir)))) {
    dir.create(ffPlot(paste0(plotSubdir)))
}

PCsToPlot = c(paste0("PC", 1:7))
for (i in seq_along(PCsToPlot)) {
    for (j in seq_along(PCsToPlot)) {
        multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                               plotCols = c(PCsToPlot[i], PCsToPlot[j]), 
                                               colorByCols=colorByCols)
        ggplot2::ggsave(filename=ffPlot(paste0(plotSubdir, "/multiColorPCAPlots_varNormPCA_", 
                                        PCsToPlot[i], "x", PCsToPlot[j], 
                                        ".pdf")), plot = multiColorPCAPlots, device = "pdf",
                        limitsize=FALSE)
        
    }
    
}


