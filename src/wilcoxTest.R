
setwd(dir = "processed/brca_PCA/analysis/sheets/old_sheets/raw_loading_average/")
a = read.csv(file = "PC_Enrichment_All_Shared_Cs_657.csv", header = TRUE, sep = ",")
erRows = grepl(pattern = "esr|eralpha", x = a[, "rsDescription"], ignore.case = TRUE)
erRows = erRows | grepl(pattern = "esr|eralpha", x = a[, "rsName"], ignore.case = TRUE)

wilcox.test(x = a$PC1[erRows], y = a$PC1[!erRows], alternative = "greater")
