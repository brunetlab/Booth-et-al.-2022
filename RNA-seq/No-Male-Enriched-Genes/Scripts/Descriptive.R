setwd("~/R-code/RNA-seq/No-Male-Enriched-Genes/")
library(Biobase)
library(RColorBrewer)

load("Data/eset_vst.RData")

pca <- prcomp(t(exprs(eset)), scale = TRUE)
s <- summary(pca)$importance[, 1:4]

#################################################################
## color coded by age and shapes for sex and presence of males ##
#################################################################
legend <- unique(pData(eset)[, c("Presence.of.males", "Age", "col", "bg", "Days.with.Males", "shp")])
legend <- legend[c(2,3,1,5,4,6), ]

pdf("Results/pca_vst.pdf",
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
legend(x=100, y=0,
       col = legend$col,
       pt.bg = legend$bg,
       pch = legend$shp,
       pt.cex = 1.3,
       legend = paste("day", legend$Age, "her (+", legend$Presence.of.males, "for", legend$Days.with.Males, "days)"),
       bty = "n")
dev.off()

pdf("Results/pca_vst_PC1-3.pdf",
    width = 6.0, height = 3.2)
par(mar = c(4.1, 4.1, 1.1, 15.1),
    xpd = TRUE)
plot(pca$x[, "PC1"],
     pca$x[, "PC3"],
     pch = pData(eset)$shp,
     cex = 1.3,
     col = pData(eset)$col,
     bg = pData(eset)$bg,
     las = 1,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC3 (", round(100*s[2,3], digits = 1),  "%)", sep = ""))
legend(x=100, y=0,
       col = legend$col,
       pt.bg = legend$bg,
       pch = legend$shp,
       pt.cex = 1.3,
       legend = paste("day", legend$Age, "her (+", legend$Presence.of.males, "for", legend$Days.with.Males, "days)"),
       bty = "n")
dev.off()


#######################################################################################################
## Hierarchical clustering, based on euclidean distance of the samples, clustering method: complete: ##
#######################################################################################################
pdf("Results/dendrogram_euclidean_complete_vst.pdf",
    width = 10, height = 5)
par(mar=c(12.1, 4.1, 1.1, 1.1))
dd <- as.dendrogram(
  hclust(
    dist(t(exprs(eset)), method = "euclidean")
  )
)
plot(dd)
dev.off()

#######################################################
## color code by batch and use shapes for treatments ##
#######################################################
col.batch <- character(length=ncol(eset))
col.batch[pData(eset)$library.batch == "1"] <- "darkturquoise"
col.batch[pData(eset)$library.batch == "2"] <- "darkviolet"
pData(eset)$col.batch <- col.batch
legend.batch <- unique(pData(eset)[, c("col.batch", "library.batch")])

pdf("Results/pca_vst_batch.pdf",
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
legend(90, 50,
       col = "grey",
       pt.bg = "grey",
       pch = legend$shp,
       pt.cex = 1.3,
       legend = paste("day", legend$Age, "her (+", legend$Presence.of.males, "for", legend$Days.with.Males, "days)"),
       bty = "n")
legend(90,100,
       col = legend.batch$col.batch,
       pt.bg = legend.batch$col.batch,
       pch = 21,
       pt.cex = 1.3,
       legend = paste("batch ", legend.batch$library.batch, sep=""),
       bty = "n")
dev.off()
