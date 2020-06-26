##########################################
##### Normalization of all data sets #####
##########################################

setwd("~/R-code/RNA-seq/No-Male-Enriched-Genes/")
library(Biobase)
library(DESeq2)
load("Data/eset_counts.RData")

pData(eset)$group <- factor(paste(pData(eset)$Days.with.Males, pData(eset)$Presence.of.males, pData(eset)$Age, sep = "_"))
dds <- DESeqDataSetFromMatrix(countData = exprs(eset),
                              design = ~ group,
                              colData = pData(eset))
dds.deseq <- DESeq(dds)
save(dds.deseq,
     file = "Data/dds_deseq.RData")


pdf("Data/dispersion_estimate.pdf")
plotDispEsts(dds.deseq)
dev.off()

vst <- varianceStabilizingTransformation(dds.deseq,
                                         blind = TRUE)
exprs(eset) <- assay(vst)
save(eset,
     file = "Data/eset_vst.RData")
write.table(exprs(eset),
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Data/vst.txt")

#############################################
#####  Normalization of only Day 3 Data #####
#############################################
setwd("~/R-code/RNA-seq/No-Male-Enriched-Genes/")
library(Biobase)
library(DESeq2)
load("Data/eset_counts.RData")

## Remove the day 7 and 13 samples from the eset
eset.day3 <- eset[, c(3:9, 22:28)]

## Normalize data from only the day 3 samples:
pData(eset.day3)$group <- factor(pData(eset.day3)$Presence.of.males)
dds.day3 <- DESeqDataSetFromMatrix(countData = exprs(eset.day3),
                                      design = ~ group,
                                      colData = pData(eset.day3))
dds.deseq.day3 <- DESeq(dds.day3)
save(dds.deseq.day3,
     file = "Data/dds_deseq_day3.RData")


pdf("Data/dispersion_estimate_day3.pdf")
plotDispEsts(dds.deseq.day3)
dev.off()

vst.day3 <- varianceStabilizingTransformation(dds.deseq.day3,
                                                 blind = TRUE)
exprs(eset.day3) <- assay(vst.day3)
save(eset.day3,
     file = "Data/eset_vst_day3.RData")
write.table(exprs(eset.day3),
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Data/vst_day3.txt")

#############################################
#####  Normalization of only Day 7 Data #####
#############################################
setwd("~/R-code/RNA-seq/No-Male-Enriched-Genes/")
library(Biobase)
library(DESeq2)
load("Data/eset_counts.RData")

## Remove the day 3 and 13 samples
eset.day7 <- eset[, c(10:21, 29:36)]

## Normalize data from only the day 7 samples:
pData(eset.day7)$group <- factor(paste(pData(eset.day7)$Days.with.Males, pData(eset.day7)$Presence.of.males, sep = "_"))
dds.day7 <- DESeqDataSetFromMatrix(countData = exprs(eset.day7),
                                   design = ~ group,
                                   colData = pData(eset.day7))
dds.deseq.day7 <- DESeq(dds.day7)
save(dds.deseq.day7,
     file = "Data/dds_deseq_day7.RData")


pdf("Data/dispersion_estimate_day7.pdf")
plotDispEsts(dds.deseq.day7)
dev.off()

vst.day7 <- varianceStabilizingTransformation(dds.deseq.day7,
                                              blind = TRUE)
exprs(eset.day7) <- assay(vst.day7)
save(eset.day7,
     file = "Data/eset_vst_day7.RData")
write.table(exprs(eset.day7),
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Data/vst_day7.txt")