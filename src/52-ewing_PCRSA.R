

setCacheDir(paste0(Sys.getenv("PROCESSED"), "ews_patients/RCache/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "ews_patients/analysis/plots/"))

########################################

source(paste0(Sys.getenv("CODE"), "/pcrsa_method_paper/src/load_process_regions_brca.R"))

#################################################################
simpleCache("bigSharedC")
# doing PCA of the methylation data
mData = bigSharedC$methylProp

simpleCache("allMPCA", {
    prcomp(t(mData), center = TRUE)
})
allMPCAWeights = as.data.table(allMPCA$rotation)


coordinates = bigSharedC[["coordinates"]]

# top10TSNE = Rtsne(X = top10MPCA$x[, 1:50], pca = FALSE, max_iter=5000,
#                   perplexity = 30)
# plot(top10TSNE$Y)
# plot(allMPCA$x[,c("PC1", "PC4")])
#plot(allMPCA$x[,c("PC2", "PC4")])
# plot(allMPCA$x[,c("PC3", "PC4")])
# taking out two samples that were outliers in PCs 2,3,4,5,
# in PC4 one was very high the other very low
ol1Ind = which.max(allMPCA$x[,c("PC4")]) # EWS_T127 (48)
ol2Ind = which.min(allMPCA$x[,c("PC4")]) # EWS_T126 (47)
olNames = colnames(bigSharedC$methylProp)[c(ol1Ind, ol2Ind)]

# post QC version
simpleCache("bigSharedC_pQC", {
    bigSharedC_pQC = bigSharedC
    #methylData = bigSharedC_pQC$methlyProp
    bigSharedC_pQC$methylProp = bigSharedC$methylProp[, !colnames(bigSharedC$methylProp) %in% olNames]
    bigSharedC_pQC$coverage = bigSharedC_pQC$coverage[, !colnames(bigSharedC$coverage) %in% olNames]
    bigSharedC_pQC
})


# then run PCA again
simpleCache("allMPCA2", {
    prcomp(t(mData)[-c(ol2Ind, ol1Ind), ], center = TRUE)
})
allMPCAWeights2 = as.data.table(allMPCA2$rotation)

mIQR = apply(mData[, -c(ol2Ind, ol1Ind)], 1, IQR)
simpleCache("top10MPCA", {
    prcomp(t(mData[mIQR >= quantile(mIQR, 0.9), -c(ol2Ind, ol1Ind)]), center = TRUE)
}, recreate = TRUE)
top10Coord = coordinates[mIQR >= quantile(mIQR, 0.9), ]
top10PCWeights = as.data.table(top10MPCA$rotation)




##################################################################
# run PC region set enrichment analysis
allMPCAString="allMPCA2"
top10MPCAString="top10MPCA"
rm("allMPCA")
PCsToAnnotate = paste0("PC", 1:10)

# using anno objects from the LOLA script
rsName = c("GSM2305313_MCF7_E2_peaks_hg38.bed", 
           lolaCoreRegionAnno$filename,
           roadmapRegionAnno$filename,
           motifRegionAnno$filename)
rsDescription = c("ER ChIPseq, MCF7 cell line, estradiol stimulation",
                  lolaCoreRegionAnno$description,
                  roadmapRegionAnno$description,
                  motifRegionAnno$description)

source(paste0(Sys.getenv("CODE"),"/aml_e3999/src/PCRSA_pipeline.R"))
# rows of mData are cytosines, cols are samples
rsEnrichResults = PCRSA_pipeline(mData = mData[, -c(ol2Ind, ol1Ind)], coordinates = coordinates, 
               GRList = GRList, 
               PCsToAnnotate=PCsToAnnotate,
               pcaCache=TRUE,
               allMPCACacheName=allMPCAString, 
               top10MPCACacheName=top10MPCAString,
               overwritePCACaches=FALSE,
               allMPCA=allMPCA2,
               top10MPCA=top10MPCA,
               rsName=rsName, rsDescription=rsDescription, 
               rsEnCache = TRUE,
               rsEnCacheName="rsEnrichment", rsEnTop10CacheName="rsEnrichmentTop10",
               overwriteResultsCaches = TRUE) 

rsEnrichment = rsEnrichResults[[1]]
rsEnrichmentTop10 = rsEnrichResults[[2]]

write.csv(x = rsEnrichment,
          file = dirData("analysis/sheets/PC_Enrichment_All_Shared_Cs_EWS.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(x = rsEnrichmentTop10,
          file = dirData("analysis/sheets/PC_Enrichment_Top_10%_Variable_Cs_EWS.csv"),
          quote = FALSE, row.names = FALSE)


# rsRankingIndex(rsEnrichment=rsEnrichment, PCsToAnnotate=PCsToAnnotate)


###################################

head(rsEnrichment[order(PC1,decreasing = TRUE)], 20)

simpleCache("rsEnrichment2", {
    rsEnrichment = pcRegionSetEnrichment(loadingMat=allMPCAWeights2, coordinateDT = coordinates, 
                                         GRList, 
                                         PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"), permute=FALSE)
    # rsNames = c("Estrogen_Receptor", loRegionAnno$filename[c(a549Ind, mcf7Ind)])
    rsNames = c("Estrogen_Receptor", loRegionAnno$filename[-sheff_dnaseInd], ewingSetNames)
    rsEnrichment[, rsNames:= rsNames]
    rsEnrichment
})

# check whether is enrichment is specific to this region set by
# seeing if loading values have a spike in the center of these region sets
# compared to surrounding genome 
GRList = lapply(GRList, resize, width = 10000, fix="center")
simpleCache("pcProfAllM", {
    pcProfAllM = pcEnrichmentProfile(loadingMat = allMPCAWeights, coordinateDT = coordinates,
                                     GRList=GRList, PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                                     binNum = 21)
    
    length(rsNames) == length(pcProfAllM)
    rsNames = c("Estrogen_Receptor", loRegionAnno$filename[-sheff_dnaseInd])
    names(pcProfAllM) <- rsNames 
    pcProfAllM
})


pcP = pcProfAllM
# plot(pcP$PC1, type="l")
# plot(pcP$PC2, type="l")
# plot(pcP$PC3, type="l")
# plot(pcP$PC4, type="l")
# plot(pcP$PC5, type="l")


grDevices::pdf(paste0(Sys.getenv("PLOTS"), "allMPCProfiles.pdf"))
for (i in 1:length(pcProfAllM)) {
    plot(pcP[[i]]$PC1, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC2, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC3, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC4, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC5, type="l") + title(rsNames[i])
}
dev.off()




# running the enrichment analysis

