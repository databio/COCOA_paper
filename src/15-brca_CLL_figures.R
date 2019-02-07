# library(projectInit)

# project.init(codeRoot = paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/"), dataDir = paste0(Sys.getenv("PROCESSED"), "brca_PCA/"))
source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)


#############################################################################
# PCA plots


############################################################################
# region set results distribution plot

simpleCache("rsEnrichmentRSMean_657", assignToVariable = "rsScores")

# summary figure of COCOA BRCA results: ER set relative ranking among region sets
rsScores$origInd = 1:nrow(rsScores)
rsScores = rsScores[order(rsScores$PC1, decreasing = TRUE), ]
erIndPC1 = grep(pattern = "suz12", x = rsScores$rsName, ignore.case = TRUE, value = TRUE)
erIndPC1 = erIndPC1 | grepl(pattern = "esr1|eralpha", x = rsScores$rsDescription, ignore.case = TRUE)
rsScores$rsDescription[erIndPC1]
plot(erIndPC1)
erIndPC1 = which(erIndPC1)
hist(erIndPC1, breaks = seq(0, 2000, by=200))

plotRSConcentration <- function(rsScores, scoreColName="PC1", 
                                colsToSearch = c("rsName", "rsDescription"), 
                                pattern, breaks) {
    
    rsScores = as.data.frame(rsScores)
    rsScores = rsScores[order(rsScores[, scoreColName], decreasing = TRUE), ]
    
    rsInd = rep(FALSE, nrow(rsScores))
    for (i in seq_along(colsToSearch)) {
        rsInd = rsInd | grepl(pattern = pattern, x = rsScores[, colsToSearch[i]], ignore.case = TRUE)
    }
    
    rsInd = which(rsInd)
    hist(rsInd, breaks = breaks)
    
}

plotRSConcentration(rsScores=rsScores, scoreColName = "PC5", pattern = "esr1|eralpha", breaks=seq(0, 2000, by=200))
plotRSConcentration(rsScores=rsScores, 
                    scoreColName = "PC1", 
                    pattern = "stat", 
                    breaks=seq(0, 2400, by=200))


# ggplot version of rs concentration
# 1 row per region set, column for rank in a given PC, 0/1 column for ER or not
rsScores

plotRSConcentration <- function(rsScores, scoreColName="PC1", 
                                colsToSearch = c("rsName", "rsDescription"), 
                                pattern, percent = FALSE) {
    # breaks

    rsRankInd = rsRankingIndex(rsScores=rsScores, PCsToAnnotate=scoreColName)

    
    rsInd = rep(FALSE, nrow(rsScores))
    for (i in seq_along(colsToSearch)) {
        rsInd = rsInd | grepl(pattern = pattern, x = rsScores[, colsToSearch[i]], ignore.case = TRUE)
    }
    
    rsScores$ofInterest = rsInd
    ofInterestDF = as.data.frame(rsInd[as.matrix(rsRankInd)])
    colnames(ofInterestDF) <- colnames(rsRankInd)
    ofInterestDF$rsRank = 1:nrow(ofInterestDF)
    categoryDistPlot = ggplot(ofInterestDF, aes(x=rsRank, weight=get(scoreColName))) + 
                              geom_histogram() + theme_classic()#+ facet_wrap(~get(scoreColName))
    return(categoryDistPlot)
    
}
##########################
# BRCA DNA methylation 

#############################
# BRCA ATAC

############################
# MOFA CLL
simpleCache("rsScore_Cor_CLL196", assignToVariable = "rsScores")
plotRSConcentration(rsScores, scoreColName="LF2", 
                                colsToSearch = c("rsName", "rsDescription"), 
                                pattern= "K562")
plotRSConcentration(rsScores, scoreColName="LF9", 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|era|gata3|foxa1|h3r17")
plotRSConcentration(rsScores, scoreColName="LF6", 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "GM12878|GM18526|GM12891|GM10847|K562|leukemia|leukaemia|lymphoma")

################################################################################
# meta region loading profile plots



################################################################################

