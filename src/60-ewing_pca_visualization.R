plotSubdir = "/"

grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir,"ewing_PC_plots.pdf"))
PCsToAnnotate = paste0("PC", 1:10)
for (i in seq_along(PCsToAnnotate)) {
    
    for (j in seq_along(PCsToAnnotate)) {
        
        plot(mPCA$x[, PCsToAnnotate[i]], mPCA$x[,PCsToAnnotate[j]], xlab=PCsToAnnotate[i], ylab=PCsToAnnotate[j])
    }
    
}
dev.off()


grDevices::pdf(paste0(Sys.getenv("PLOTS"), plotSubdir,"ewing_PC1_PC2.pdf"))


plot(mPCA$x[, "PC1"], mPCA$x[, "PC2"], xlab="PC1", ylab="PC2")

dev.off()
