## Here, I am comparing the DEseq2 results of the glp-1 datasets (this paper) to the N2 and fem-1 results from 
## the self-sperm paper (Booth et al. eLife 2019) using Spearman's

rm(list=ls())

library(Biobase)
library(pheatmap)
library(RColorBrewer)
##library(DESeq2)
setwd("~/R-code/RNA-seq/No-Male-Enriched-Genes/")

################################################
### Create data frames of the DEseq2 results ###
################################################
## DEGs from the male-factors paper:
N2malesvctrl.3 <- read.delim("Results/DEseq/Groups/N2malesvctrl.3.txt", 
                             stringsAsFactors = FALSE,
                             row.names = 1) 
N2malesvctrl.7.1 <- read.delim("Results/DEseq/Groups/N2malesvctrl.7.1.txt", 
                               stringsAsFactors = FALSE,
                               row.names = 1) 
N2malesvctrl.7.5 <- read.delim("Results/DEseq/Groups/N2malesvctrl.7.5.txt", 
                               stringsAsFactors = FALSE,
                               row.names = 1)

identical(rownames(N2malesvctrl.3), rownames(N2malesvctrl.7.1))
identical(rownames(N2malesvctrl.3), rownames(N2malesvctrl.7.5))
MF.results <- data.frame(FC.N2malesvctrl.3 = N2malesvctrl.3$log2FoldChange,
                         padj.N2malesvctrl.3 = N2malesvctrl.3$padj,
                         FC.N2malesvctrl.7.1 = N2malesvctrl.7.1$log2FoldChange,
                         padj.N2malesvctrl.7.1 = N2malesvctrl.7.1$padj,
                         FC.N2malesvctrl.7.5 = N2malesvctrl.7.5$log2FoldChange,
                         padj.N2malesvctrl.7.5 = N2malesvctrl.7.5$padj,
                         gene = rownames(N2malesvctrl.7.1))
rownames(MF.results) <- rownames(N2malesvctrl.7.1)

## DEGs from the Self-sperm paper:
N3MvN3U <- read.delim("Data/2018-05-MID-N2fem-1/DEseq/N3MvN3U.txt",
                      stringsAsFactors = FALSE,
                      row.names = 1)
F3MvF3U <- read.delim("Data/2018-05-MID-N2fem-1/DEseq/F3MvF3U.txt",
                      stringsAsFactors = FALSE,
                      row.names = 1)

identical(rownames(F3MvF3U), rownames(N3MvN3U))
SS.results <- data.frame(FC.F3MvF3U = F3MvF3U$log2FoldChange,
                         padj.F3MvF3U = F3MvF3U$padj,
                         FC.N3MvN3U = N3MvN3U$log2FoldChange,
                         padj.N3MvN3U = N3MvN3U$padj,
                         gene = rownames(F3MvF3U))
rownames(SS.results) <- rownames(F3MvF3U)

## Merge the two datasets:
all.results <- merge(MF.results, SS.results, by.x = "gene", by.y = "gene") 
rownames(all.results) <- all.results$gene
comps <- c("N2malesvctrl.3", "N2malesvctrl.7.1", "N2malesvctrl.7.5", "N3MvN3U", "F3MvF3U")

## Remove rows without a padj value (na):
all.results <- all.results[!(is.na(all.results$padj.N2malesvctrl.3)),]
all.results <- all.results[!(is.na(all.results$padj.N2malesvctrl.7.1)),]
all.results <- all.results[!(is.na(all.results$padj.N2malesvctrl.7.5)),]
all.results <- all.results[!(is.na(all.results$padj.N3MvN3U)),]
all.results <- all.results[!(is.na(all.results$padj.F3MvF3U)),]
## 8291 genes

## Remove pseudogenes
load("Data/eset_counts.RData")
PGs <- subset(eset, fData(eset)$biotype == "pseudogene")
pseudogenes <- fData(PGs)$gene.name
all.results <- all.results[!is.element(rownames(all.results), pseudogenes), ]
## 8182 genes after removing pseudogenes

## Remove the "gene" column and the columns of the padj values
all.results <- all.results[ -c(1) ]
all.results.FC <- all.results[ -c(2,4,6,8,10)]



## Create a matrix of Pearsons's correlations of the fold-changes
corr.pearson.matrix <- cor(all.results.FC, method = c("pearson"))

## Plot this as a heatmap
breaklist = seq(0, 1, length.out = 100)

pdf("Results/N2-fem-1-comparison/Pearsons.pdf",
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
