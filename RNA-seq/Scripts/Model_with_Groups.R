## Simplest model: Just define age/presence.of.males groups: GE ~ groups
## HEADS UP: If p values are reported as 0, that just means that they are actually p < 1/num.samples and should be reported like this.
## in the mean time...it is accurate to say that p<0.001


setwd("~/R-code/RNA-seq/")
library(Biobase)
library(DESeq2)
load("Data/eset_counts.RData")

pData(eset)$group <- factor(paste(pData(eset)$Sex, pData(eset)$Age, pData(eset)$Presence.of.males, pData(eset)$Days.with.Males, sep = "_"))
dds <- DESeqDataSetFromMatrix(countData = exprs(eset),
                             design = ~ group,
                             colData = pData(eset))
dds.deseq <- DESeq(dds)


## The numbers of differentially expressed genes for key comparisons:
## control 7v3:
sum(results(dds.deseq, contrast = list("groupHer_7_none_0", "groupHer_3_none_0"))$padj < 0.05, na.rm = TRUE) ## 4363
## control 13v7:
sum(results(dds.deseq, contrast = list("groupHer_13_none_0", "groupHer_7_none_0"))$padj < 0.05, na.rm = TRUE) ## 3091
## +N2 males, age and length of exposure:
sum(results(dds.deseq, contrast = list("groupHer_7_N2.males_5", "groupHer_3_N2.males_1"))$padj < 0.05, na.rm = TRUE) ##  6509
## +N2 males, age only:
sum(results(dds.deseq, contrast = list("groupHer_7_N2.males_1", "groupHer_3_N2.males_1"))$padj < 0.05, na.rm = TRUE)
## N2 vs control at day 3:
sum(results(dds.deseq, contrast = list("groupHer_3_N2.males_1", "groupHer_3_none_0"))$padj < 0.05, na.rm = TRUE)
## N2 vs control at day 7, longer exposure:
sum(results(dds.deseq, contrast = list("groupHer_7_N2.males_5", "groupHer_7_none_0"))$padj < 0.05, na.rm = TRUE)
## N2 vs control at day 7, short exposure:
sum(results(dds.deseq, contrast = list("groupHer_7_N2.males_1", "groupHer_7_none_0"))$padj < 0.05, na.rm = TRUE)
## +N2 at day 7, long v. short exposure:
sum(results(dds.deseq, contrast = list("groupHer_7_N2.males_5", "groupHer_7_N2.males_1"))$padj < 0.05, na.rm = TRUE)
## N2 males versus day 3 glp-1 her:
sum(results(dds.deseq, contrast = list("groupMale_7_N2.males_5", "groupHer_3_none_0"))$padj < 0.05, na.rm = TRUE)

## Saving the results of the DEseq2 analysis:
N2males7.5v3.1 <- results(dds.deseq, contrast = list("groupHer_7_N2.males_5", "groupHer_3_N2.males_1"))
N2males7.1v3.1 <- results(dds.deseq, contrast = list("groupHer_7_N2.males_1", "groupHer_3_N2.males_1"))
ctrl7v3 <- results(dds.deseq, contrast = list("groupHer_7_none_0", "groupHer_3_none_0"))
save(N2males7.5v3.1, N2males7.1v3.1, ctrl7v3,
     file = "Data/7v3_deseq.RData")

ctrl13v3 <- results(dds.deseq, contrast = list("groupHer_13_none_0", "groupHer_3_none_0"))
ctrl13v7 <- results(dds.deseq, contrast = list("groupHer_13_none_0", "groupHer_7_none_0"))
save(ctrl7v3, ctrl13v3, ctrl13v7,
     file = "Data/Her_only_aging_deseq.RData")

N2malesvctrl.3 <- results(dds.deseq, contrast = list("groupHer_3_N2.males_1", "groupHer_3_none_0"))
N2malesvctrl.7.5 <- results(dds.deseq, contrast = list("groupHer_7_N2.males_5", "groupHer_7_none_0"))
N2malesvctrl.7.1 <- results(dds.deseq, contrast = list("groupHer_7_N2.males_1", "groupHer_7_none_0"))
save(N2malesvctrl.3, N2malesvctrl.7.5, N2malesvctrl.7.1, 
     file = "Data/males_v_nomales_deseq.RData")

N2males7.5v7.1 <- results(dds.deseq, contrast = list("groupHer_7_N2.males_5", "groupHer_7_N2.males_1"))
save(N2males7.5v7.1,
     file = "Data/long_v_short_males7_deseq.RData")

N2malevher3 <- results(dds.deseq, contrast = list("groupMale_7_N2.males_5", "groupHer_3_none_0"))
N2malevher7 <- results(dds.deseq, contrast = list("groupMale_7_N2.males_5", "groupHer_7_none_0"))
N2malevher13 <- results(dds.deseq, contrast = list("groupMale_7_N2.males_5", "groupHer_13_none_0"))
save(N2malevher3, N2malevher7, N2malevher13,
     file = "Data/male_v_her_deseq.RData")

all.comparisons <- list(N2males7.5v3.1 = N2males7.5v3.1,
                        N2males7.1v3.1 = N2males7.1v3.1,
                        ctrl7v3 = ctrl7v3,
                        ctrl13v3 = ctrl13v3,
                        ctrl13v7 = ctrl13v7,
                        N2malesvctrl.3 = N2malesvctrl.3, 
                        N2malesvctrl.7.5 = N2malesvctrl.7.5, 
                        N2malesvctrl.7.1 = N2malesvctrl.7.1,
                        N2males7.5v7.1 = N2males7.5v7.1,
                        N2malevher3 = N2malevher3, 
                        N2malevher7 = N2malevher7, 
                        N2malevher13 = N2malevher13)

for(i in 1:length(all.comparisons)){
  comparison.i <- all.comparisons[[i]]
  file.name <- names(all.comparisons[i])
  file.path <- paste0("Results/DEseq/Groups/", file.name, ".txt")
  write.table(comparison.i,
              col.names = NA,
              sep = "\t",
              quote = FALSE,
              file = file.path)
}

