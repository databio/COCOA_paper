# fit distribution to permutations to approximate p value with 
# small number of permutations

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
devtools::load_all(ffCode("COCOA/"))

nCores = 1 # detectCores() - 1
options("mc.cores"=nCores)

scriptID = "21-testPermDist"
plotSubdir = "21-testPermDist/"
dataID = "CLL196MOFA_cov"
sheetsDir = ffProc("COCOA_paper/analysis/sheets/")

if (!dir.exists(ffPlot(plotSubdir))) {
    dir.create(ffPlot(plotSubdir))
}

set.seed(1234)
nPerm = 300

#############################################################################

load(paste0(getCacheDir(),"rsPermScores_", dataID, ".RData"))
simpleCache(paste0("rsScore_", dataID), assignToVariable = "realRSScores")
########
nullDistList = lapply(X = 1:nrow(rsPermScores[[1]]),
                      FUN = function(x) extractNullDist(resultsList=rsPermScores, rsInd = x))
nullDistList = lapply(nullDistList, FUN=function(x) as.data.frame(x)[, colsToAnnotate])

for (j in c(1, 58, 1100, 1500, 2000)) {
    rsInd = j
    a = fitGammaNullDistr(nullDistList[[rsInd]], method = "mle", force = TRUE)
    pdf(ffPlot(paste0(plotSubdir, "gammaMLEFitExamples", rsInd, ".pdf")))
    for (i in seq_along(colsToAnnotate)) {
        plot(a[[colsToAnnotate[i]]])    
    }
    
    dev.off()
    
    b = fitGammaNullDistr(nullDistList[[rsInd]], method="mme")
    pdf(ffPlot(paste0(plotSubdir, "gammaMMEFitExamples", rsInd, ".pdf")))
    for (i in seq_along(colsToAnnotate)) {
        plot(b[[colsToAnnotate[i]]])    
    }
    
    dev.off()
}


#####
gPValDF = getGammaPVal(scores = realRSScores[, colsToAnnotate], nullDistList = nullDistList)
gPValDF = cbind(gPValDF, realRSScores[, colnames(realRSScores)[!(colnames(realRSScores) %in% colsToAnnotate)]])

multiNiceHist(file = ffPlot(paste0(plotSubdir, "pValDist", dataID, ".pdf")), dataDF = -log10(gPValDF[colsToAnnotate]),
              colsToPlot = colsToAnnotate, xLabels = "p-value",
              binwidth = 1, boundary = 0, yLabel = "Number of region sets",
              plotTitles = paste0("Distribution of region set p-values (-log10), ", colsToAnnotate),
              ggExpr = "+ylim(0, 2270)")
