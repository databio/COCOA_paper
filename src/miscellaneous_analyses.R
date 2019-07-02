# analyses/scripts that do not need to be published with the paper

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


