# create survival plots

library(survival)


dataID = "kircMethyl"

#############################################################################
simpleCache(paste0("rsScores_", dataID, "Cor"), assignToVariable = "realRSScores")

topRSInd = rsRankingIndex(rsScores = realRSScores, signalCol = "cancerStage")$cancerStage


# look at average correlation per region to determine whether all CpGs 
# are going in same (important for aggregating since a simple average
# of the DNA methylation will run into problems if some CpGs go in opposite directions)
averagePerRegion(signal = , signalCoord = , regionSet = GRList[[topRSInd[5]]], absVal = FALSE)

# first test whether DNA methylation is associated with cancer stage 


# test whether DNA methylation is associated with survival
# kaplan meier plot: groups are predicted good outcome vs predicted bad outcome
patSurv = Surv(patientMetadata_pqc$os_months, rep(1, length(patientMetadata_pqc$os_months)))
kmFit = survfit(patSurv ~ 1)
patSurv2 = Surv(patientMetadata_pqc$os_months, rep_len(c(1,1,0), length.out = length(patientMetadata_pqc$os_months)))
kmFit2 = survfit(formula = patSurv2 ~ patientMetadata_pqc$Complex, data = patientMetadata_pqc)

plot(kmFit2)


###### cox proportional hazards model

# covariates
covariateData = patientMetadata_pqc

# add ML model predictions to covariate DF
# join/merge?

patSurv = Surv(patientMetadata_pqc$os_months, rep(1, length(patientMetadata_pqc$os_months)))
coxph(patSurv ~ age + model_pred + NPM1, data = covariateData)