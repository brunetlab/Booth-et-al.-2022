## I am comparing my RNA-seq DEseq2 results and Cheng's glp-1 her x WT male microarray results
## using Venn diagrams and hypergeometric test

rm(list=ls())
library(Vennerable)
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
                             padj.N2malesvctrl.3 = N2malesvctrl.3$padj,
                             FC.N2malesvctrl.7.1 = N2malesvctrl.7.1$log2FoldChange,
                             padj.N2malesvctrl.7.1 = N2malesvctrl.7.1$padj,
                             FC.N2malesvctrl.7.5 = N2malesvctrl.7.5$log2FoldChange,
                             padj.N2malesvctrl.7.5 = N2malesvctrl.7.5$padj,
                             gene = rownames(N2malesvctrl.7.1))
rownames(RNAseq.results) <- rownames(N2malesvctrl.7.1)

## DEGs from the Cheng's microarrays:
microarray <- read.delim("Data/WTmale-microarray-FC.txt",
                         stringsAsFactors = FALSE,
                         row.names = 1)

microarray.results <- data.frame(gene = rownames(microarray),
                                 FC.microarray = microarray$average.FC)
rownames(microarray.results) <- rownames(microarray)

## Cheng's SAM results do not have the calculated p-values for each gene
## Instead, only the genes that meet the FDR cut-off are listed (all are significant)
## As a work-around for the script, I will include these genes as significant by giving them an arbitraty adjusted p-value of 0.01 
## this way, they can be included in this analysis pipeline

## 1. create a list of the significant genes using the SAM results:
microarray.sig <- read.delim("Data/Sig.Microarray.Results.WTmales.txt",
                             stringsAsFactors = FALSE,
                             row.names = NULL)
## 2. Include significance in the microarray.results dataframe by adding a column of arbitrary significant (0.01) and n.s. (1) values
microarray.results$padj.microarray[is.element(rownames(microarray.results), microarray.sig[[1]])] <- "0.01"
microarray.results$padj.microarray[!is.element(rownames(microarray.results), microarray.sig[[1]])] <- "1"


## Merge the two datasets (only keeping genes that are found in both the microarray and RNA-seq):
all.results <- merge(RNAseq.results, microarray.results, by.x = "gene", by.y = "gene") 
rownames(all.results) <- all.results$gene
## Remove the "gene" column
all.results <- all.results[ -c(1) ]
comps <- c("N2malesvctrl.3", "N2malesvctrl.7.1", "N2malesvctrl.7.5", "microarray")

## Identify genes that are significant for each dataset:
p <- 0.05
fc <- 0
## positive:
sig.results.positive <- list()
for(i in 1:length(comps)){
  comp.i <- comps[i]
  p.col <- paste0("padj.", comp.i)
  FC.col <- paste0("FC.", comp.i)
  ind <- which(all.results[, p.col] < p & all.results[, FC.col] > fc )
  sig.results.positive[[i]] <- rownames(all.results)[ind]
}
names(sig.results.positive) <- comps

## negative:
sig.results.negative <- list()
for(i in 1:length(comps)){
  comp.i <- comps[i]
  p.col <- paste0("padj.", comp.i)
  FC.col <- paste0("FC.", comp.i)
  ind <- which(all.results[, p.col] < p & all.results[, FC.col] < fc)
  sig.results.negative[[i]] <- rownames(all.results)[ind]
}
names(sig.results.negative) <- comps

## Create a Venn Diagram
## only the 7.5 RNA-seq condition and the microarray comparison are in the figures (Supp. Fig. 3b)
Venn.results.pos <- Venn(sig.results.positive)
Venn.results.pos.array.3.1 <- Venn.results.pos[ , c("N2malesvctrl.3",
                                           "microarray")]
Venn.results.pos.array.7.1 <- Venn.results.pos[ , c("N2malesvctrl.7.1",
                                            "microarray")]

Venn.results.pos.array.7.5 <- Venn.results.pos[ , c("N2malesvctrl.7.5",
                                            "microarray")]

Venn.results.neg <- Venn(sig.results.negative)
Venn.results.neg.array.3.1 <- Venn.results.neg[ , c("N2malesvctrl.3",
                                                    "microarray")]
Venn.results.neg.array.7.1 <- Venn.results.neg[ , c("N2malesvctrl.7.1",
                                                    "microarray")]

Venn.results.neg.array.7.5 <- Venn.results.neg[ , c("N2malesvctrl.7.5",
                                                    "microarray")]

pdf("Venn.results.pos.array-3.1.pdf",
    width = 7, height = 7)
plot(Venn.results.pos.array.3.1,
     doWeights = TRUE,
     type = "circles")
dev.off()

pdf("Venn.results.pos.array-7.1.pdf",
    width = 7, height = 7)
plot(Venn.results.pos.array.7.1,
     doWeights = TRUE,
     type = "circles")
dev.off()

pdf("Venn.results.pos.array-7.5.pdf",
    width = 7, height = 7)
plot(Venn.results.pos.array.7.5,
     doWeights = TRUE,
     type = "circles")
dev.off()

pdf("Venn.results.neg.array-3.1.pdf",
    width = 7, height = 7)
plot(Venn.results.neg.array.3.1,
     doWeights = TRUE,
     type = "circles")
dev.off()

pdf("Venn.results.neg.array-7.1.pdf",
    width = 7, height = 7)
plot(Venn.results.neg.array.7.1,
     doWeights = TRUE,
     type = "circles")
dev.off()

pdf("Venn.results.neg.array-7.5.pdf",
    width = 7, height = 7)
plot(Venn.results.neg.array.7.5,
     doWeights = TRUE,
     type = "circles")
dev.off()


## Hypergeometric test
## N is the total population size (# of genes detected in all datasets)
## M is the number of successes in the population (significant genes in the microarray)
## k is the sample size (# of significant genes in the RNA-seq dataset)
## q is the number of successes in the sample (# of genes significant in both the microarray and the RNA-seq)
N <- length(rownames(all.results))
M.pos <- length(sig.results.positive[["microarray"]])
M.neg <- length(sig.results.negative[["microarray"]])
N.pos <- N - M.pos
N.neg <- N - M.neg


hyper.positive <- list()
for(i in 1:length(comps)){
  comp.i <- comps[i]
  k <- length(sig.results.positive[[comp.i]])
  q <- length(intersect(sig.results.positive[[comp.i]], sig.results.positive[["microarray"]]))
  hyper.positive[[i]] <- phyper(q-1, M.pos, N.pos, k, lower.tail = FALSE, log.p = FALSE)
}
names(hyper.positive) <- comps
write.csv(hyper.positive, file = "posDEGs-hypergeometric.csv")

hyper.negative <- list()
for(i in 1:length(comps)){
  comp.i <- comps[i]
  k <- length(sig.results.negative[[comp.i]])
  q <- length(intersect(sig.results.negative[[comp.i]], sig.results.negative[["microarray"]]))
  hyper.negative[[i]] <- phyper(q-1, M.neg, N.neg, k, lower.tail = FALSE, log.p = FALSE)
}
names(hyper.negative) <- comps
write.csv(hyper.negative, file = "negDEGs-hypergeometric.csv")
