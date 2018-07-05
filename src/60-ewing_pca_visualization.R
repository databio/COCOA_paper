plotSubdir = "/"
library(ggplot2)

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

pcDF = data.frame(PC1=mPCA$x[, "PC1"], PC2 = mPCA$x[, "PC2"])

pc1_2 = ggplot(data = pcDF, mapping = aes(x=PC1, y=PC2)) + geom_point(size=3) + theme_classic()
ggsave(paste0(Sys.getenv("PLOTS"), plotSubdir,"ewing_PC_plots.pdf"))
