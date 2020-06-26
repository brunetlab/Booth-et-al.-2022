## I am comparing my RNA-seq DEseq2 results and Cheng's glp-1 her x WT male microarray results using Pearson's
## This results is presented in Supplemental Figure 3a

rm(list=ls())

library(Biobase)
library(pheatmap)
library(RColorBrewer)
##library(DESeq2)
setwd("~/R-code/Microarrays/")

###############################################
### Create data frame of the DEseq2 results ###
###############################################
## DEGs from my RNA-seq results:
## Using the DEseq2 results with removal of the male-enriched genes 
N2malesvctrl.3 <- read.delim("Data/RNA-seq/N2malesvctrl.3.txt", 
                             stringsAsFactors = FALSE,
                             row.names = 1) 
N2malesvctrl.7.1 <- read.delim("Data/RNA-seq/N2malesvctrl.7.1.txt", 
                               stringsAsFactors = FALSE,
                               row.names = 1) 
N2malesvctrl.7.5 <- read.delim("Data/RNA-seq/N2malesvctrl.7.5.txt", 
                               stringsAsFactors = FALSE,
                               row.names = 1)

identical(rownames(N2malesvctrl.3), rownames(N2malesvctrl.7.1))
identical(rownames(N2malesvctrl.3), rownames(N2malesvctrl.7.5))
RNAseq.results <- data.frame(FC.N2malesvctrl.3 = N2malesvctrl.3$log2FoldChange,
                         FC.N2malesvctrl.7.1 = N2malesvctrl.7.1$log2FoldChange,
                         FC.N2malesvctrl.7.5 = N2malesvctrl.7.5$log2FoldChange,
                         gene = rownames(N2malesvctrl.7.1))
rownames(RNAseq.results) <- rownames(N2malesvctrl.7.1)

## DEGs from the Cheng's microarrays:
microarray <- read.delim("Data/WTmale-microarray-FC.txt",
                      stringsAsFactors = FALSE,
                      row.names = 1)

microarray.results <- data.frame(gene = rownames(microarray),
                                 FC.microarray = microarray$average.FC)
rownames(microarray.results) <- rownames(microarray)

## Merge the two datasets:
all.results <- merge(RNAseq.results, microarray.results, by.x = "gene", by.y = "gene") 
rownames(all.results) <- all.results$gene
## Remove the "gene" column
all.results <- all.results[ -c(1) ]

## Create a matrix of Pearsons's correlations of the fold-changes
corr.pearson.matrix <- cor(all.results, method = c("pearson"))

## Plot this as a heatmap
breaklist = seq(0, 1, length.out = 100)

pdf("RNAseq-microarray-Pearsons.pdf",
    width = 4, height = 3.5, onefile=FALSE)
pheatmap(corr.pearson.matrix,
         scale = "none",
         color = colorRampPalette(brewer.pal(n = 5, name = "Purples"))(length(breaklist)),
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         border_color = NA,
         breaks = breaklist)
dev.off()
