rm(list=ls())

setwd("~/R-code/Microarrays/")

## Here I will compare the array results from Cheng to each other
## These comparisons should reveal the role of male sperm vs. seminal fluid
## I am using the resutls from Cheng's 7 glp-1 her x WT males (day 2 and 3) and 6 glp-1 her x fer-6 males SAM results

## Results from Cheng
arrays.results <- read.delim("WTmales_fer-6males.txt",
                             stringsAsFactors = FALSE,
                             row.names = NULL)

sig.results <- list(glp1hxfer6mUP = arrays.results[["glp1hxfer6mUP"]],
                    glp1hxfer6mDOWN = arrays.results[["glp1hxfer6mDOWN"]], 
                    glp1hxWTmUP = arrays.results[["glp1hxWTmUP"]], 
                    glp1hxWTmDOWN = arrays.results[["glp1hxWTmDOWN"]])

## Creates a file with the genes listed in the comparisons (sperm-dependent):
common.UP <- intersect(sig.results[["glp1hxfer6mUP"]], sig.results[["glp1hxWTmUP"]])
write(common.UP,
      file = "Lists/common.WTm-fer6m-UP.txt")

common.DOWN <- intersect(sig.results[["glp1hxfer6mDOWN"]], sig.results[["glp1hxWTmDOWN"]])
write(common.DOWN,
      file = "Lists/common.WTm-fer6m-DOWN.txt")

WTmalespecific.UP <- setdiff(sig.results[["glp1hxWTmUP"]], sig.results[["glp1hxfer6mUP"]])
write(WTmalespecific.UP,
      file = "Lists/WTmalespecific_UP.txt")

WTmalespecific.DOWN <- setdiff(sig.results[["glp1hxWTmDOWN"]], sig.results[["glp1hxfer6mDOWN"]])
write(WTmalespecific.DOWN,
      file = "Lists/WTmalespecific_DOWN.txt")

fer6malespecific.UP <- setdiff(sig.results[["glp1hxfer6mUP"]], sig.results[["glp1hxWTmUP"]])
write(fer6malespecific.UP,
      file = "Lists/fer6malespecific_UP.txt")

fer6malespecific.DOWN <- setdiff(sig.results[["glp1hxfer6mDOWN"]], sig.results[["glp1hxWTmDOWN"]])
write(fer6malespecific.DOWN,
      file = "Lists/fer6malespecific_DOWN.txt")

common <- c(common.UP, common.DOWN)
WTonly <- c(WTmalespecific.UP, WTmalespecific.DOWN)


###########################################
## Create heat maps with the subset data ##
###########################################
library(Biobase)
library(pheatmap)
library(RColorBrewer)

## Set color scaling 
breaklist = seq(-3, 3, length.out = 100)

## fold-changes for each microarray replicate from Cheng's data:
arrays.FC <- read.delim("microarrayFC.txt",
                             stringsAsFactors = FALSE,
                             row.names = NULL)
rownames(arrays.FC) <- arrays.FC$gene.name

## subset the genes changing in hermaphrodites with both WT males and fer-6 males (seminal fluid induced)
SF.genes <- subset(arrays.FC, arrays.FC$gene.ID %in% common)
dim(SF.genes)
## 130 15
## Remove first two columns (have gene names)
SF.genes <- subset(SF.genes, select = -c(1,2))
dim(SF.genes)
## 130 13
## Remove the WT males results to cluster and plot only the fer6 results:
SF.genes <- subset(SF.genes, select = -c(1:7))

##heatmap of seminal fluid genes (common genes changing in hermaphrodites in presence WT males and fer-6 males )
pdf("Heatmaps/SFgenes-common-WTm_fer6m.pdf",
    width = 5, height = 10, onefile=FALSE)
pheatmap(SF.genes,
         scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaklist)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         border_color = NA,
         breaks = breaklist)
dev.off()


## subset the genes changing in hermaphrodites with only WT males and not with fer-6 males (sperm induced?)
sperm.genes <- subset(arrays.FC, arrays.FC$gene.ID %in% WTonly)
dim(sperm.genes)
## 1341 15
## Remove first two columns (have gene names)
sperm.genes <- subset(sperm.genes, select = -c(1,2))
dim(sperm.genes)
## 1341 13
## Remove columns with the fer-6 males to plot only the WT male results
sperm.genes <- subset(sperm.genes, select = c(1:7))

##heatmap of sperm genes (common genes changing in hermaphrodites in presence WT males but not fer-6 males )
pdf("Heatmaps/spermgenes-only-WTm.pdf",
    width = 5, height = 10, onefile=FALSE)
pheatmap(sperm.genes,
         scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaklist)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         border_color = NA,
         breaks = breaklist)
dev.off()


## Create a heatmap of the SAM-identified significant male-conditioned media hermaphrodite genes
## fold-changes for each microarray replicate from Cheng's data:
MCParrays.FC <- read.delim("MCPmicroarrayFC.txt",
                        stringsAsFactors = FALSE,
                        row.names = NULL)
rownames(MCParrays.FC) <- MCParrays.FC$gene.name
MCParrays.FC <- subset(MCParrays.FC, select = -c(1,2))

pdf("Heatmaps/MCPgenes.pdf",
    width = 5, height = 10, onefile=FALSE)
pheatmap(MCParrays.FC,
         scale = "none",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaklist)),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         border_color = NA,
         breaks = breaklist)
dev.off()

## Create a Venn diagram of the microarray results:
library(Vennerable)

## Add MCP gene names listed by up or down-regulated to the sig.results lists:
MCP <- read.delim("MCP-male_v_her.txt",
                  stringsAsFactors = FALSE,
                  row.names = 1)

sig.results[["MCP.UP"]] <- MCP$her.UP
sig.results[["MCP.DOWN"]] <- MCP$her.DOWN
## Remove blank elements from the list:
sig.results <- lapply(sig.results, function(x) x[x != ""])

Venn.results <- Venn(sig.results)
Venn.UP <- Venn.results[ , c("glp1hxWTmUP", "glp1hxfer6mUP", "MCP.UP")]
Venn.DOWN <- Venn.results[ , c("glp1hxWTmDOWN", "glp1hxfer6mDOWN", "MCP.DOWN")]

pdf("Venn.microarrays.UP.pdf",
    width = 7, height = 7)
plot(Venn.UP,
     doWeights = TRUE,
     type = "circles")
dev.off()

pdf("Venn.microarrays.DOWN.pdf",
    width = 7, height = 7)
plot(Venn.DOWN,
     doWeights = TRUE,
     type = "circles")
dev.off()



