rm(list=ls())

library(Biobase)
library(pheatmap)
library(RColorBrewer)
##library(DESeq2)
setwd("~/R-code/RNA-seq/No-Male-Enriched-Genes/Results/")

####################################################################################
### Create data frame to identify sets of genes based on differential expression ###
####################################################################################
N2malesvctrl.3 <- read.delim("DEseq/Groups/N2malesvctrl.3.txt", 
                             stringsAsFactors = FALSE,
                             row.names = 1) 
N2malesvctrl.7.1 <- read.delim("DEseq/Groups/N2malesvctrl.7.1.txt", 
                               stringsAsFactors = FALSE,
                               row.names = 1) 
N2malesvctrl.7.5 <- read.delim("DEseq/Groups/N2malesvctrl.7.5.txt", 
                               stringsAsFactors = FALSE,
                               row.names = 1)

identical(rownames(N2malesvctrl.3), rownames(N2malesvctrl.7.1))
identical(rownames(N2malesvctrl.3), rownames(N2malesvctrl.7.5))
all.results <- data.frame(FC.N2malesvctrl.3 = N2malesvctrl.3$log2FoldChange,
                          padj.N2malesvctrl.3 = N2malesvctrl.3$padj,
                          FC.N2malesvctrl.7.1 = N2malesvctrl.7.1$log2FoldChange,
                          padj.N2malesvctrl.7.1 = N2malesvctrl.7.1$padj,
                          FC.N2malesvctrl.7.5 = N2malesvctrl.7.5$log2FoldChange,
                          padj.N2malesvctrl.7.5 = N2malesvctrl.7.5$padj)
rownames(all.results) <- rownames(N2malesvctrl.7.1)
comps <- c("N2malesvctrl.3", "N2malesvctrl.7.1", "N2malesvctrl.7.5")

## Remove pseudogenes
load("../Data/eset_counts.RData")
PGs <- subset(eset, fData(eset)$biotype == "pseudogene")
pseudogenes <- fData(PGs)$gene.name
all.results <- all.results[!is.element(rownames(all.results), pseudogenes), ]

## Create a list of genes that are significant in each group (p<0.05)
p <- 0.05
## Significant in Either Direction:
sig.results <- list()
for(i in 1:length(comps)){
  comp.i <- comps[i]
  p.col <- paste0("padj.", comp.i)
  FC.col <- paste0("FC.", comp.i)
  ind <- which(all.results[, p.col] < p)
  sig.results[[i]] <- rownames(all.results)[ind]
}
names(sig.results) <- comps

##########################################################################################
## Create data frame of only the fold change values to be used to input to the heatmaps ##
##########################################################################################
FC.results <- data.frame(N2malesvctrl.3 = N2malesvctrl.3$log2FoldChange,
                         N2malesvctrl.7.1 = N2malesvctrl.7.1$log2FoldChange,
                         N2malesvctrl.7.5 = N2malesvctrl.7.5$log2FoldChange)
rownames(FC.results) <- rownames(N2malesvctrl.7.1)

## Remove pseudogenes:
FC.results <- FC.results[!is.element(rownames(FC.results), pseudogenes), ]


#########################
### Set color scaling ###
#########################
breaklist = seq(-5, 5, length.out = 100)



##########################################################################################
### Only genes significantly changing in at least two time points in response to males ###
##########################################################################################
## Significant genes are from DEseq (see above)
atleast2.common <- intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.1"]])
atleast2.common <- append(atleast2.common, intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.5"]]))
atleast2.common <- append(atleast2.common, intersect(sig.results[["N2malesvctrl.7.1"]], sig.results[["N2malesvctrl.7.5"]]))
atleast2.common <- unique(atleast2.common)

## Heatmap:
atleast2.common.FC <- subset(FC.results, rownames(FC.results) %in% atleast2.common)
dim(atleast2.common.FC)
## 81  3

pdf("DEseq/Groups/Heatmaps/heatmap_log2FC_3.1-7.1-7.5_2-way-intersection.pdf",
    width = 3, height = 10, onefile=FALSE)
pheatmap(atleast2.common.FC,
         scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaklist)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         border_color = NA,
         breaks = breaklist)
dev.off()


###########################################################################
### Only genes significantly changing in both ages in response to males ###
###########################################################################
## Significant genes are from DEseq (see above)
age.common <- intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.1"]])
age.common <- append(age.common, intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.5"]]))
age.common <- unique(age.common)

## Heatmap:
genes.common.FC <- subset(FC.results, rownames(FC.results) %in% age.common)
dim(genes.common.FC)
## 10  3

pdf("DEseq/Groups/Heatmaps/heatmap_log2FC_3.1-7.1-7.5_age-independent.pdf",
    width = 3, height = 3.5, onefile=FALSE)
pheatmap(genes.common.FC,
         scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaklist)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         border_color = NA,
         breaks = breaklist)
dev.off()
