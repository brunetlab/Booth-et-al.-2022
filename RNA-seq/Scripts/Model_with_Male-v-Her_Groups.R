## Simplest model: Just define age/presence.of.males groups: GE ~ groups
## HEADS UP: If p values are reported as 0, that just means that they are actually p < 1/num.samples and should be reported like this.
## in the mean time...it is accurate to say that p<0.001


setwd("~/R-code/RNA-seq/")
library(Biobase)
library(DESeq2)
load("Data/eset_counts.RData")

## Subset data to include only the male samples and the day 3 and 7 hermaphrodites with out males
eset.sex <- eset[, c(3:21)]

## I will combine the day 3 and 7 samples with each other
pData(eset.sex)$sex <- factor(paste(pData(eset.sex)$Sex))
dds <- DESeqDataSetFromMatrix(countData = exprs(eset.sex),
                             design = ~ sex,
                             colData = pData(eset.sex))
dds.deseq <- DESeq(dds)


## The numbers of differentially expressed genes for key comparisons:
## males versus hermaphrodites:
sum(results(dds.deseq, contrast = list("sexMale", "sexHer"))$padj < 0.05, na.rm = TRUE) ## 10224
sum(results(dds.deseq, contrast = list("sexMale", "sexHer"))$padj < 0.01, na.rm = TRUE) ## 8841


## Saving the results of the DEseq2 analysis:
male.versus.herm3.7 <- results(dds.deseq, contrast = list("sexMale", "sexHer"))
save(male.versus.herm3.7, file = "Data/male_v_her37_deseq.RData")
write.table(male.versus.herm3.7,
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Results/DEseq/Male-v-Her_Groups/male_v_her37.txt")

## Saving just the male enriched ones:
ind <- which(male.versus.herm3.7$padj < 0.05 & male.versus.herm3.7$log2FoldChange > 0 )
sig.results.positive <- rownames(male.versus.herm3.7)[ind]
save(sig.results.positive, file = "Data/male-enriched-p05.RData")

ind.str <- which(male.versus.herm3.7$padj < 0.01 & male.versus.herm3.7$log2FoldChange > 0 )
sig.results.positive.stringent <- rownames(male.versus.herm3.7)[ind.str]
save(sig.results.positive.stringent, file = "Data/male-enriched-p01.RData")

