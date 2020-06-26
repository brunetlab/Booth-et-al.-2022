setwd("~/R-code/RNA-seq/No-Male-Enriched-Genes/")
library(Biobase)
library(RColorBrewer)

load("Data/eset_vst_day7.RData")
eset <- eset.day7

pca <- prcomp(t(exprs(eset)), scale = FALSE)
s <- summary(pca)$importance[, 1:4]

#################################################################
## color coded by age and shapes for sex and presence of males ##
#################################################################
legend <- unique(pData(eset)[, c("Presence.of.males", "col", "bg", "shp", "Days.with.Males")])

pdf("Results/pca_vst_day7.pdf",
    width = 6.0, height = 3.2)
par(mar = c(4.1, 4.1, 1.1, 15.1),
    xpd = TRUE)
plot(pca$x[, "PC1"],
     pca$x[, "PC2"],
     pch = pData(eset)$shp,
     cex = 1.3,
     col = pData(eset)$col,
     bg = pData(eset)$bg,
     las = 1,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""))
legend(x=70, y=0,
       col = legend$col,
       pt.bg = legend$bg,
       pch = legend$shp,
       pt.cex = 1.3,
       legend = paste("day 7 her (+", legend$Presence.of.males, "for", legend$Days.with.Males, "day)"),
       bty = "n")
dev.off()


#######################################################
## color code by batch and use shapes for treatments ##
#######################################################
col.batch <- character(length=ncol(eset))
col.batch[pData(eset)$library.batch == "1"] <- "darkturquoise"
col.batch[pData(eset)$library.batch == "2"] <- "darkviolet"
pData(eset)$col.batch <- col.batch
legend.batch <- unique(pData(eset)[, c("col.batch", "library.batch")])

pdf("Results/pca_vst_batch_day7.pdf",
    width = 5.9, height = 3.2)
par(mar = c(4.1, 4.1, 1.1, 14.6),
    xpd = TRUE)
plot(pca$x[, "PC1"],
     pca$x[, "PC2"],
     pch = pData(eset)$shp,
     cex = 1.3,
     col = pData(eset)$col.batch,
     bg = pData(eset)$col.batch,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""))
legend(100, 0,
       col = "grey",
       pt.bg = "grey",
       pch = legend$shp,
       pt.cex = 1.3,
       legend = paste("day 7 her (+", legend$Presence.of.males, "for", legend$Days.with.Males, "day)"),
       bty = "n")
legend(100,20,
       col = legend.batch$col.batch,
       pt.bg = legend.batch$col.batch,
       pch = 21,
       pt.cex = 1.3,
       legend = paste("batch ", legend.batch$library.batch, sep=""),
       bty = "n")
dev.off()
