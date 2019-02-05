
library(caret)

source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))

setwd(paste0(Sys.getenv("PROCESSED"), "COCOA_paper/analysis/"))
patientMetadata = brcaMetadata # already screened out patients with incomplete ER or PGR mutation status
# there should be 657 such patients
set.seed(1234)


########################################################################

simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList")
brcaMetadata = as.data.frame(brcaMetadata)
row.names(brcaMetadata) = brcaMetadata$subject_ID

# patientMetadata should have already screened out patients without ER/PGR status
# resulting in 657 patients
hasER_PGR_IDs = patientMetadata[, subject_ID]
filteredMData = brcaMList[["methylProp"]][, 
                                          colnames(brcaMList[["methylProp"]]) %in% hasER_PGR_IDs] 

methylData = t(filteredMData)
erStatus = brcaMetadata[row.names(methylData), "ER_status"]
erStatus = as.factor(erStatus)
methylData = cbind(methylData, erStatus)



trainingParams = trainControl(method = "none", summaryFunction = twoClassSummary)

# random forest
simpleCache("erStatusRF", {
   erStatusRF = train(erStatus ~ ., data=methylData, 
                    method = "ranger", trControl = trainingParams)
   erStatusRF
})

