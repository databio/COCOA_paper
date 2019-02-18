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



# plotRSConcentration <- function(rsScores, scoreColName="PC1", 
#                                 colsToSearch = c("rsName", "rsDescription"), 
#                                 pattern, breaks) {
#     
#     rsScores = as.data.frame(rsScores)
#     rsScores = rsScores[order(rsScores[, scoreColName], decreasing = TRUE), ]
#     
#     rsInd = rep(FALSE, nrow(rsScores))
#     for (i in seq_along(colsToSearch)) {
#         rsInd = rsInd | grepl(pattern = pattern, x = rsScores[, colsToSearch[i]], ignore.case = TRUE)
#     }
#     
#     rsInd = which(rsInd)
#     hist(rsInd, breaks = breaks)
#     
# }
# 
# plotRSConcentration(rsScores=rsScores, scoreColName = "PC5", pattern = "esr1|eralpha", breaks=seq(0, 2000, by=200))
# plotRSConcentration(rsScores=rsScores, 
#                     scoreColName = "PC1", 
#                     pattern = "stat", 
#                     breaks=seq(0, 2400, by=200))



##########################
# BRCA DNA methylation 
simpleCache("rsEnrichmentRSMean_657", assignToVariable = "rsScores")

# summary figure of COCOA BRCA results: ER set relative ranking among region sets
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|eralpha") + ggtitle("Estrogen receptor region sets")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|eralpha|gata3|foxa1|h3r17") + ggtitle("Estrogen receptor-related region sets")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|eralpha|gata3|foxa1|h3r17")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "ezh2|suz12")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "pol2")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k9")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k4me3")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k36")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "k562")
plotRSConcentration(rsScores, scoreColName=paste0("PC", 1:9), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "mcf7|mcf-7")
#############################
# BRCA ATAC

############################
# MOFA CLL
simpleCache("rsScore_Cor_CLL196", assignToVariable = "rsScores")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                                colsToSearch = c("rsName", "rsDescription"), 
                                pattern= "K562")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|era|gata3|foxa1|h3r17")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "GM12878|GM18526|GM12891|GM10847|K562|leukemia|leukaemia|lymphoma")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k9")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k4me1")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k4me3")
plotRSConcentration(rsScores, scoreColName=c(paste0("LF", c(1:3, 5:7, 9))), 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "h3k36")
################################################################################
# meta region loading profile plots



################################################################################

