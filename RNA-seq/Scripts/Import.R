## Make expressionSet with count matrix.

setwd("~/R-code/RNA-seq/")
source("Scripts/make.transparent.R")
library(Biobase)

files <- list.files("Data/Star-ReadCounts",
                    pattern = "ReadsPerGene",
                    full.names = TRUE)
samples <- gsub("Data/Star-ReadCounts/", "", files)
samples <- gsub("_ReadsPerGene.out.tab", "", samples)
system(paste("wc -l", files[1]))
## 46752 Data/Star-ReadCounts/LB-1-TAAGGCGA-TAAGGCGA_S1_ReadsPerGene.out.tab
count.matrix <- matrix(nrow = 46752,
                       ncol = length(files),
                       dimnames = list(NULL, samples))

## counts for unstranded RNA-seq in 2nd column.
for(i in seq(along = files)){
  data.i <- read.table(files[i],
                       colClasses = c("NULL", "integer","NULL", "NULL"),
                       sep = "\t")[,1]
  count.matrix[, samples[i]] <- data.i
  rm(data.i)
}
genes <- read.table(files[i],
                    colClasses = c("character", "NULL","NULL", "NULL"),
                    sep = "\t")[,1]
rownames(count.matrix) <- genes

## Get Sample info:
pData <- read.table("Data/SampleInfo.txt",
                    sep = "\t",
                    row.names = 1,
                    header = TRUE)
pData <- pData[match(colnames(count.matrix), rownames(pData)), ]
pData$col[pData$Presence.of.males == "none"
          & pData$Age == 3] <- "deepskyblue"
pData$col[pData$Presence.of.males == "none"
          & pData$Age == 7] <- "firebrick1"
pData$col[pData$Presence.of.males == "none"
          & pData$Age == 13] <- "darkgoldenrod2"
pData$col[pData$Presence.of.males == "N2 males"
          & pData$Age == 3] <- "dodgerblue3"
pData$col[pData$Presence.of.males == "N2 males"
          & pData$Age == 7] <- "firebrick"
pData$col[pData$Sex == "Male"] <- "darkviolet"
pData$bg <- make.transparent(pData$col, alpha = 200)
pData$shp[pData$Presence.of.males == "none"] <- 21
pData$shp[pData$Presence.of.males == "N2 males"
          & pData$Days.with.Males == "1"
          & pData$Sex == "Her"] <- 24
pData$shp[pData$Presence.of.males == "N2 males"
          & pData$Days.with.Males == "5"
          & pData$Sex == "Her"] <- 23
pData$shp[pData$Sex == "Male"] <- 10
pData$Presence.of.males <- relevel(pData$Presence.of.males, ref = "none")
pData$fileName <- rownames(pData)


## Get gene info:
tab <- read.table("Data/gene_info.txt",
                  sep = ";",
                  stringsAsFactors = FALSE)
fData <- data.frame(gene.id = c(rownames(count.matrix)[1:4], gsub("gene_id ", "", tab$V1)),
                    gene.name = c(rownames(count.matrix)[1:4], gsub(" gene_name ", "", tab$V2)),
                    source = c(rownames(count.matrix)[1:4], gsub(" gene_source ", "", tab$V3)),
                    biotype = c(rownames(count.matrix)[1:4], gsub(" gene_biotype ", "", tab$V4)),
                    stringsAsFactors = FALSE)
rownames(fData) <- fData$gene.name
identical(fData$gene.id, rownames(count.matrix))
## TRUE

eset <- ExpressionSet(assayData = count.matrix)
pData(eset) <- pData
fData(eset) <- fData

rownames(eset) <- fData(eset)$gene.name
colnames(eset) <- pData(eset)$Sample.Name

eset <- eset[, order(colnames(eset))]

## Remove the first 4 rows (unmapped reads etc):
eset <- eset[-c(1:4),]

## Filter out genes with low coverage across samples.
## I will only keep genes that have a cpm (counts per million) value of at least 1 in at least 3 samples:
cov.per.sample <- colSums(exprs(eset))
norm.fact <- rep(cov.per.sample, each = nrow(eset))
cpm <- exprs(eset) / norm.fact * 10^6
ind.keep <- rowSums(cpm >= 1) >= 3
table(ind.keep)
## FALSE  TRUE
## 30042 16706
genes.keep <- names(ind.keep)[ind.keep]
eset <- eset[genes.keep, ]

table(fData(eset)$biotype)
##        antisense        lincRNA          miRNA          ncRNA          piRNA protein_coding     pseudogene 
##            13             73              7            141             20          16013           402
##        rRNA         snoRNA          snRNA           tRNA 
##         5            28              3              1 
dim(eset)
## Features  Samples
##     16706       40

## Coverage per sample:
quantile(colSums(exprs(eset)))
##       0%      25%      50%      75%     100%
##  3650980  9797989 12360667 18272383 31137429

## Number of genes per sample:
quantile(colSums(exprs(eset) != 0))
##      0%     25%       50%     75%    100%
##  14330.00 14945.75 15493.50 16243.75 16551.00

save(eset,
     file = "Data/eset_counts.RData")

