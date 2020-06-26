## This script will create a heatmap to summarize ChIP Atlas results for my screen hit genes
## Only selected TFs are presented in figure 3
## The results from each gene-TF combo are independent of each other

rm(list=ls())
library(Biobase)
library(pheatmap)
library(RColorBrewer)
setwd("~/R-code/Microarrays/")

chip.results <- read.delim("Data/Screen-genes-chip-atlas.txt", 
                             stringsAsFactors = FALSE,
                             row.names = 1)


pdf("Screen-genes-chip-MACS2-scores.pdf",
    width = 8, height = 5, onefile=FALSE)
pheatmap(chip.results,
         scale = "none",
         color = colorRampPalette(brewer.pal(n = 5, name = "RdPu"))(50),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames = TRUE,
         na_col = "grey95",
         border_color = NA)
dev.off()