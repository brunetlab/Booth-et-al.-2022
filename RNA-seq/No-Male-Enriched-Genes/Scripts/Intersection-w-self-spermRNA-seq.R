rm(list=ls())

library(Biobase)
library(pheatmap)
library(RColorBrewer)
##library(DESeq2)
setwd("~/R-code/RNA-seq/No-Male-Enriched-Genes/")

#####################################################################################
### Create data frames to identify sets of genes based on differential expression ###
#####################################################################################
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
## 10,352 genes in MF.results

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
## 14,147 genes in SS.results

## Merge the two datasets:
all.results <- merge(MF.results, SS.results, by.x = "gene", by.y = "gene") 
rownames(all.results) <- all.results$gene
comps <- c("N2malesvctrl.3", "N2malesvctrl.7.1", "N2malesvctrl.7.5", "N3MvN3U", "F3MvF3U")
## 8936 genes after merge

## Remove pseudogenes
load("Data/eset_counts.RData")
PGs <- subset(eset, fData(eset)$biotype == "pseudogene")
pseudogenes <- fData(PGs)$gene.name
all.results <- all.results[!is.element(rownames(all.results), pseudogenes), ]
## 8813 genes after removing pseudogenes

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
FC.results <- data.frame(N2malesvctrl.3 = all.results$FC.N2malesvctrl.3,
                         N2malesvctrl.7.1 = all.results$FC.N2malesvctrl.7.1,
                         N2malesvctrl.7.5 = all.results$FC.N2malesvctrl.7.5,
                         N3MvN3U = all.results$FC.N3MvN3U,
                         F3MvF3U = all.results$FC.F3MvF3U)
rownames(FC.results) <- rownames(all.results)

#########################
### Set color scaling ###
#########################
breaklist = seq(-5, 5, length.out = 100)

#########################################################
### Only genes significantly changing in all datasets ###
#########################################################
all.common <- intersect(sig.results[["F3MvF3U"]], 
                        intersect(sig.results[["N3MvN3U"]], 
                                  intersect(sig.results[["N2malesvctrl.7.5"]], 
                                            intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.1"]]))))
write.table(all.common, sep = "\t", quote = FALSE, file = "Results/N2-fem-1-comparison/all-common-genes.txt")

## Heatmap:
all.common.FC <- subset(FC.results, rownames(FC.results) %in% all.common)
pdf("Results/N2-fem-1-comparison/heatmap_log2FC_allcommongenes.pdf",
    width = 3, height = 3, onefile=FALSE)
pheatmap(all.common.FC,
         scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaklist)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         border_color = NA,
         breaks = breaklist)
dev.off()


################################################################
### Only genes significantly changing in at least 4 datasets ###
################################################################
four.common <- intersect(sig.results[["N3MvN3U"]],
                         intersect(sig.results[["N2malesvctrl.7.5"]],
                                   intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.1"]])))
four.common <- append(four.common, intersect(sig.results[["F3MvF3U"]],
                                             intersect(sig.results[["N2malesvctrl.7.5"]],
                                                       intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.1"]]))))
four.common <- append(four.common, intersect(sig.results[["F3MvF3U"]],
                                             intersect(sig.results[["N3MvN3U"]],
                                                       intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.1"]]))))
four.common <- append(four.common, intersect(sig.results[["F3MvF3U"]],
                                             intersect(sig.results[["N3MvN3U"]],
                                                       intersect(sig.results[["N2malesvctrl.7.5"]], sig.results[["N2malesvctrl.7.1"]]))))
four.common <- append(four.common, intersect(sig.results[["F3MvF3U"]],
                                             intersect(sig.results[["N3MvN3U"]],
                                                       intersect(sig.results[["N2malesvctrl.7.5"]], sig.results[["N2malesvctrl.3"]]))))
four.common <- unique(four.common)
write.table(four.common, sep = "\t", quote = FALSE, file = "Results/N2-fem-1-comparison/four-way-common-genes.txt")

## Heatmap:
four.common.FC <- subset(FC.results, rownames(FC.results) %in% four.common)
pdf("Results/N2-fem-1-comparison/heatmap_log2FC_four-way-commongenes.pdf",
    width = 3, height = 6, onefile=FALSE)
pheatmap(four.common.FC,
         scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaklist)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         border_color = NA,
         breaks = breaklist)
dev.off()

################################################################
### Only genes significantly changing in at least 3 datasets ###
################################################################
three.common <- intersect(sig.results[["N2malesvctrl.7.5"]], intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.1"]]))
three.common <- append(three.common, intersect(sig.results[["N3MvN3U"]], intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.1"]])))
three.common <- append(three.common, intersect(sig.results[["F3MvF3U"]], intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.1"]])))
three.common <- append(three.common, intersect(sig.results[["N3MvN3U"]], intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.5"]])))
three.common <- append(three.common, intersect(sig.results[["F3MvF3U"]], intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.5"]])))
three.common <- append(three.common, intersect(sig.results[["N3MvN3U"]], intersect(sig.results[["N2malesvctrl.3"]], sig.results[["F3MvF3U"]])))
three.common <- append(three.common, intersect(sig.results[["N3MvN3U"]], intersect(sig.results[["N2malesvctrl.7.1"]], sig.results[["N2malesvctrl.7.5"]])))
three.common <- append(three.common, intersect(sig.results[["N3MvN3U"]], intersect(sig.results[["N2malesvctrl.7.1"]], sig.results[["F3MvF3U"]])))
three.common <- append(three.common, intersect(sig.results[["F3MvF3U"]], intersect(sig.results[["N2malesvctrl.7.1"]], sig.results[["N2malesvctrl.7.5"]])))
three.common <- append(three.common, intersect(sig.results[["N3MvN3U"]], intersect(sig.results[["F3MvF3U"]], sig.results[["N2malesvctrl.7.5"]])))
three.common <- unique(three.common)
write.table(three.common, sep = "\t", quote = FALSE, file = "Results/N2-fem-1-comparison/three-way-common-genes.txt")

## Heatmap:
three.common.FC <- subset(FC.results, rownames(FC.results) %in% three.common)
pdf("Results/N2-fem-1-comparison/heatmap_log2FC_three-way-commongenes.pdf",
    width = 3, height = 20, onefile=FALSE)
pheatmap(three.common.FC,
         scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaklist)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         border_color = NA,
         breaks = breaklist)
dev.off()


################################################################
### Only genes significantly changing in at least 2 datasets ###
################################################################
two.common <- intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.1"]])
two.common <- append(two.common, intersect(sig.results[["N2malesvctrl.3"]], sig.results[["N2malesvctrl.7.5"]]))
two.common <- append(two.common, intersect(sig.results[["N2malesvctrl.7.1"]], sig.results[["N2malesvctrl.7.5"]]))
two.common <- append(two.common, intersect(sig.results[["N3MvN3U"]], sig.results[["F3MvF3U"]]))
two.common <- append(two.common, intersect(sig.results[["N3MvN3U"]], sig.results[["N2malesvctrl.3"]]))
two.common <- append(two.common, intersect(sig.results[["N3MvN3U"]], sig.results[["N2malesvctrl.7.1"]]))
two.common <- append(two.common, intersect(sig.results[["N3MvN3U"]], sig.results[["N2malesvctrl.7.5"]]))
two.common <- append(two.common, intersect(sig.results[["F3MvF3U"]], sig.results[["N2malesvctrl.3"]]))
two.common <- append(two.common, intersect(sig.results[["F3MvF3U"]], sig.results[["N2malesvctrl.7.1"]]))
two.common <- append(two.common, intersect(sig.results[["F3MvF3U"]], sig.results[["N2malesvctrl.7.5"]]))

two.common <- unique(two.common)
write.table(two.common, sep = "\t", quote = FALSE, file = "Results/N2-fem-1-comparison/two-way-common-genes.txt")

## Heatmap:
two.common.FC <- subset(FC.results, rownames(FC.results) %in% two.common)
pdf("Results/N2-fem-1-comparison/heatmap_log2FC_two-way-commongenes.pdf",
    width = 3, height = 85, onefile=FALSE)
pheatmap(two.common.FC,
         scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaklist)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         border_color = NA,
         breaks = breaklist)
dev.off()


