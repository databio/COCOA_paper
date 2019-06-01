# load libraries and prepare environment for other scripts

library(LOLA)
library(simpleCache)
library(data.table)
library(ggplot2)
library(GenomicRanges) # GRangesList, resize
# library(caret)
library(RGenomeUtils)
library(MIRA)
library(ComplexHeatmap)
library(gridExtra) #marrangeGrob for colorClusterPlots()
# some of the environmental variables from aml/.../00-init.R will need to be reset
source(paste0(Sys.getenv("CODE"), "aml_e3999/src/00-genericFunctions.R" )) 
library(MultiAssayExperiment)
library(folderfun)
# the AML init script will set a different plots directory
# Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/plots/"))

# source(paste0(Sys.getenv("CODE"), "COCOA/R/COCOA.R"))
# source(paste0(Sys.getenv("CODE"), "COCOA/R/visualization.R"))
library(COCOA)

brcaMetadata = fread(paste0(Sys.getenv("CODE"), 
                         "COCOA_paper/metadata/brca_metadata.csv"))
# only keep patients who have definitive status for ER and PGR
brcaMetadata = brcaMetadata[brcaMetadata$ER_status %in% 
                                                     c("Positive", "Negative"), ]
brcaMetadata = brcaMetadata[brcaMetadata$PGR_status %in% 
                                c("Positive", "Negative"), ]


# indexed clinical data (brca_clinical_metadata.tsv) has more up to date follow up info
brcaMetadata2 = read.table(file = paste0(Sys.getenv("CODE"), 
                                         "COCOA_paper/metadata/brca_clinical_metadata.tsv"), 
                           sep = "\t", header = TRUE)

patientMetadata = merge(brcaMetadata, brcaMetadata2[, c("bcr_patient_barcode", 
                                                       "vital_status", 
                                                       "days_to_death", 
                                                       "days_to_last_follow_up")],
                        by.x="subject_ID", 
                        by.y="bcr_patient_barcode", all.x=TRUE)

row.names(patientMetadata) <- patientMetadata$subject_ID

process_brca_expr = function(exprDF) {
    
}


setSubPlotDir = function(dirName) {
    Sys.setenv("SUBPLOTDIR"=paste0(dirName, "/"))
}

dirPlot = function(plotFile) {
    if (!is.na(Sys.getenv("SUBPLOTDIR", unset = NA))) {
        return(paste0(Sys.getenv("PLOTS"), Sys.getenv("SUBPLOTDIR"), plotFile))    
    } else(
        return(paste0(Sys.getenv("PLOTS"), plotFile))
    )
    
}
setff("Plot", Sys.getenv("PLOTS"))

setff("Code", paste0(Sys.getenv("CODE")))
setff("ProjCode", paste0(Sys.getenv("CODE"), "COCOA_paper/"))
dirCode = function(.file="") {
    return(paste0(Sys.getenv("CODE"), "COCOA_paper/", .file))
}



# set environment
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/plots/"))
setCacheDir(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/RCache/"))

############ functions to add to COCOA ##############################
# ggplot version of rs concentration
# 1 row per region set, column for rank in a given PC, 0/1 column for ER or not

# TODO: add check that scoreColNames are present as columns of rsScores
plotRSConcentration <- function(rsScores, scoreColName="PC1", 
                                colsToSearch = c("rsName", "rsDescription"), 
                                pattern, percent = FALSE, 
                                binwidth=50) {
    # breaks
    rsScores = as.data.frame(rsScores)
    rsRankInd = rsRankingIndex(rsScores=rsScores, PCsToAnnotate=scoreColName)
    
    
    rsInd = rep(FALSE, nrow(rsScores))
    for (i in seq_along(colsToSearch)) {
        rsInd = rsInd | grepl(pattern = pattern, x = rsScores[, colsToSearch[i]], ignore.case = TRUE)
    }
    
    rsScores$ofInterest = rsInd
    # reorder logical vector for each PC
    ofInterestDF = as.data.frame(mapply(FUN = function(ind) rsInd[as.matrix(rsRankInd)[, ind]], 
                          ind=seq_along(scoreColName))) 
    colnames(ofInterestDF) <- colnames(rsRankInd)
    ofInterestDF$rsRank = 1:nrow(ofInterestDF)
    # reshaping to long format
    ofInterestDF = tidyr::gather(ofInterestDF, key="PC", value="of_interest", scoreColName)
    categoryDistPlot = ggplot(ofInterestDF, aes(x=rsRank, weight=of_interest)) + 
        geom_histogram(binwidth=binwidth, boundary=0, closed="right") + theme_classic() + xlab("Region set rank") +
        ylab(paste0("Number of region sets (binwidth=", binwidth, ")")) + facet_wrap(~PC)
    return(categoryDistPlot)
}


# @param dataMat columns of dataMat should be samples/patients, rows should be genomic signal
# (each row corresponds to one genomic coordinate/range)
# @param featureMat Rows should be samples, columns should be "features" 
# (whatever you want to get correlation with: eg PC scores),
# all columns in featureMat will be used (subset when passing to function
# in order to not use all columns)
# @param center logical object. Should rows in dataMat be centered based on
# their means? (subtracting row means from each row)
#
# If a row in dataMat has 0 stand. deviation, correlation will be set to 0
# instead of NA as would be done by cor()
#
# returns a matrix where rows are the genomic signal (eg a CpG or region) and
# columns are the columns of featureMat
createCorFeatureMat = function(dataMat, featureMat, 
                               centerDataMat=TRUE, centerFeatureMat = TRUE) {
    
    featureMat = as.matrix(featureMat)
    
    if (centerDataMat) {
        cpgMeans = rowMeans(dataMat)
        # centering before calculating correlation
        dataMat = apply(X = dataMat, MARGIN = 2, function(x) x - cpgMeans)
        
    }
    
    if (centerFeatureMat) {
        featureMeans = colMeans(featureMat)
        # centering before calculating correlation
        featureMat = t(apply(X = featureMat, MARGIN = 1, function(x) x - featureMeans))
        featureMat = as.matrix(featureMat)
    }
    
    
    # create feature correlation matrix with PCs (rows: features/CpGs, columns:PCs)
    # how much do features correlate with each PC?
    featurePCCor = apply(X = featureMat, MARGIN = 2, function(y) apply(X = dataMat, 1, FUN = function(x) cor(x = x, y)))
    
    # if standard deviation of the data was zero, NA will be produced
    # set to 0 because no standard deviation means no correlation with attribute of interest
    featurePCCor[is.na(featurePCCor)] = 0
    colnames(featurePCCor) <- colnames(featureMat)
    
    return(featurePCCor)
    # corLoadRatio = loadingMat[, PCsToAnnotate] / featurePCCor 
    # hist(corLoadRatio[, "PC10"])
}

################## not sure whether to add this one to COCOA
# convenience function to quickly make meta-region loading profiles
# first calculates profiles, then normalizes and plots
# output is marrangeGrob and should be saved with ggsave
makeMetaRegionPlots <- function(loadingMat, signalCoord, GRList, rsNames, PCsToAnnotate, binNum) {
    
    pcProf = lapply(X = GRList, function(x) getLoadingProfile(loadingMat = loadingMat, 
                                                              signalCoord = signalCoord, 
                                                              regionSet = x, PCsToAnnotate = PCsToAnnotate,
                                                              binNum = binNum))
    
    pcP = copy(pcProf)
    pcP = lapply(pcP,FUN = as.data.table)
    
    
    # average loading value from each PC to normalize so PCs can be compared with each other
    avLoad = apply(X = loadingMat[, PCsToAnnotate], MARGIN = 2, FUN = function(x) mean(abs(x)))
    
    # normalize
    # pcP = lapply(pcP, FUN = function(x) t(apply(X = x, MARGIN = 1, FUN = function(y) y - c(0, avLoad))))
    pcP = lapply(pcP, FUN = function(x) x[, mapply(FUN = function(y, z) get(y) - z, y=PCsToAnnotate, z = avLoad)])
    pcP = lapply(pcP, FUN = function(x) data.table(regionGroupID=1:nrow(x), x))
    
    # for the plot scale
    maxVal = max(sapply(pcP, FUN = function(x) max(x[, .SD, .SDcols=PCsToAnnotate])))
    minVal = min(sapply(pcP, FUN = function(x) min(x[, .SD, .SDcols=PCsToAnnotate])))
    
    # convert to long format for plots
    pcP = lapply(X = pcP, FUN = function(x) tidyr::gather(data = x, key = "PC", value="loading_value", PCsToAnnotate))
    pcP = lapply(X = pcP, as.data.table)
    pcP = lapply(pcP, function(x) x[, PC := factor(PC, levels = PCsToAnnotate)])
    
    # stack overflow for wrapping plot title
    wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n") 
    
    
    profilePList = list()
    for (i in seq_along(pcP)) {
        
        thisRS = pcP[[i]]
        
        
        
        profilePList[[i]] = ggplot(data = thisRS, mapping = aes(x =regionGroupID , y = loading_value)) + 
            geom_line() + ylim(c(minVal, maxVal)) + facet_wrap(facets = "PC") + 
            ggtitle(label = wrapper(rsNames[i], width=30)) + xlab("Genome around Region Set, 14 kb") + 
            ylab("Normalized Loading Value") + 
            theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
        profilePList[[i]]
        # plot(pcP[[i]]$PC1, type="l") + title(rsNames[i])
        # 
        
        #xLabels = xAxisForRegionPlots2()
        
    }
    multiProfileP = marrangeGrob(profilePList, ncol = 2, nrow = 2)
    
    return(multiProfileP)
}


########## MOFA/CLL analysis #########################
# gets coordinates for CLL methyl
# returns a list with methylProp that has methylation and "methylCoord"
# that has corresponding genomic coordinates
prepareCLLMethyl = function() {
    
    # get microarray data
    eh = ExperimentHub()
    meth = eh[[names(query(eh, "CLLmethylation"))]] # EH1071
    # rows are cpgs, columns are samples
    methData = assay(meth)
    dataProbeNames = row.names(methData)
    
    
    # match probe names with genomic coordinates
    # ls('package:FDb.InfiniumMethylation.hg19')
    m450kAnno = get450k()
    length(m450kAnno)
    
    # get coordinates in same order as CLL data
    methCoord = m450kAnno[dataProbeNames]
    all(names(methCoord) == dataProbeNames)
    methCoordDT = COCOA:::grToDt(methCoord)
    # keep start coordinate as CpG site
    methCoordDT = methCoordDT[, .(chr, start)]
    methCoord = COCOA:::dtToGr(methCoordDT)
    if (nrow(methCoordDT) != nrow(methData)) {
        stop("error matching probes to coordinates")
    }
    
    return(list(methylProp = methData, methylCoord = methCoord))
    
}