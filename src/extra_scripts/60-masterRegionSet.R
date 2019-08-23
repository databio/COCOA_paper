# combining all region sets into a master list of regions

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))
###########################################################
# reading in the region sets
# load LOLA database
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))

#################################################################
allMPCAString = "allMPCA_657" #  "allMPCA_657"

# simpleCache(allMPCAString, assignToVariable = "mPCA")  
# brcaCoord = mPCA
simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList")
brcaCoord = brcaMList$coordinates

allRegions = do.call("c", GRList)
length(allRegions) # 243041837
allRegionsMerged = reduce(x = allRegions)
length(allRegionsMerged) # 42566
median(width(allRegionsMerged)) # 316
hist(log(width(allRegionsMerged), base = 10))
allRegionsMerged[which.max((width(allRegionsMerged)))]
allRegionsMerged = allRegionsMerged[order(width(allRegionsMerged), decreasing = TRUE)]
head(allRegionsMerged)
sum(width(allRegionsMerged)) # 2,913,545,469

# removing very large regions
allRegionsMerged = reduce(x = allRegions[width(allRegions) <= 1e4])
length(allRegionsMerged) # 60600
median(width(allRegionsMerged)) # 945
hist(log(width(allRegionsMerged), base = 10))
width(allRegionsMerged[which.max((width(allRegionsMerged)))])
allRegionsMerged = allRegionsMerged[order(width(allRegionsMerged), decreasing = TRUE)]
head(width(allRegionsMerged))
sum(width(allRegionsMerged)) # 2,889,101,946


# segmenting the genome into regions
3e9 / 500

# how many much genome regions that overlap with DNA methylation data
brcaRegions= getOLRegions(GRangesList(allRegions), COCOA:::dtToGr(brcaCoord))[[1]]
length(brcaRegions) # 32202915
hist(log(width(brcaRegions), base = 10))
brcaRegionsMerged = reduce(brcaRegions)
length(brcaRegionsMerged) # 6493
median(width(brcaRegionsMerged)) # 23519
hist(log(width(brcaRegionsMerged), base = 10))
brcaRegionsMerged[which.max((width(brcaRegionsMerged)))]
brcaRegionsMerged = brcaRegionsMerged[order(width(brcaRegionsMerged), decreasing = TRUE)]
head(brcaRegionsMerged)
sum(width(brcaRegionsMerged)) # 2,554,986,098
2554986098 / 3e9

# getting rid of regions that are very large (and perhaps less informative)
hist(log(width(brcaRegions), base = 10))
sum(width(brcaRegions) > 1e4)
sum((width(brcaRegions) <= 1e4) & (width(brcaRegions) > 5000))

# screen out regions larger than 1e4
brcaRegionsSmaller = brcaRegions[width(brcaRegions) <= 1e4]

hist(log(width(brcaRegionsSmaller), base = 10))
brcaRegionsMergedSmaller = reduce(brcaRegionsSmaller)
length(brcaRegionsMergedSmaller) # 58000
median(width(brcaRegionsMergedSmaller)) # 8624
hist(log(width(brcaRegionsMergedSmaller), base = 10))
width(brcaRegionsMergedSmaller[which.max((width(brcaRegionsMergedSmaller)))])
brcaRegionsMergedSmaller = brcaRegionsMergedSmaller[order(width(brcaRegionsMergedSmaller), decreasing = TRUE)]
head(brcaRegionsMergedSmaller)
brcaSCov = sum(width(brcaRegionsMergedSmaller)) # 622,443,703
brcaSCov / 3e9
brcaSCov / 500 # 1,244,887
3e9 / 500

#################################################################################
# similar analysis with Sheffield DNAse data

lolaPath1 = paste0(Sys.getenv("RESOURCES"), "/regions/LOLACore/hg38/")
regionSetDB = loadRegionDB(lolaPath1)
loRegionAnno = regionSetDB$regionAnno
lolaCoreRegionAnno = loRegionAnno
# a549Ind = grep("a549", loRegionAnno$cellType, ignore.case = TRUE)
sheff_dnaseInd = grep("sheffield_dnase", loRegionAnno$collection, ignore.case = TRUE)

lolaCoreRegionAnno = lolaCoreRegionAnno[sheff_dnaseInd, ]
GRList1 = GRangesList(regionSetDB$regionGRL[sheff_dnaseInd])


allRegions = do.call("c", GRList1)
length(allRegions) # 2888031
hist(log10(width(allRegions)))
max(width(allRegions)) # 2,450,903
allRegionsMerged = reduce(x = allRegions)
length(allRegionsMerged) # 2,885,104
median(width(allRegionsMerged)) # 150
hist(log(width(allRegionsMerged), base = 10))
width(allRegionsMerged[which.max((width(allRegionsMerged)))])
allRegionsMerged = allRegionsMerged[order(width(allRegionsMerged), decreasing = TRUE)]
head(width(allRegionsMerged))
sheffDNAseCov = sum(width(allRegionsMerged)) # 435,966,549
sheffDNAseCov / 500 # 871933.1

sum(width(allRegions) > 10e4)
allRegionsMerged = reduce(x = allRegions[width(allRegions) <= 10e4])
length(allRegionsMerged) # 2,886,387
median(width(allRegionsMerged)) # 150
hist(log(width(allRegionsMerged), base = 10))
width(allRegionsMerged[which.max((width(allRegionsMerged)))]) # 22913
allRegionsMerged = allRegionsMerged[order(width(allRegionsMerged), decreasing = TRUE)]
head(width(allRegionsMerged))
sheffDNAseCov = sum(width(allRegionsMerged)) # 433,573,414
sheffDNAseCov / 500 # 867146.8



