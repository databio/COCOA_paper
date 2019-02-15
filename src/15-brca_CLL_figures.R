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
plotRSConcentration(rsScores, scoreColName="PC1", 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|eralpha|gata3|foxa1|h3r17")
plotRSConcentration(rsScores, scoreColName="PC4", 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|eralpha|gata3|foxa1|h3r17")
plotRSConcentration(rsScores, scoreColName="PC4", 
                    colsToSearch = c("rsName", "rsDescription"), 
                    pattern= "esr|eralpha|gata3|foxa1|h3r17")


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

