## Make expressionSet with count matrix.
## module load destiny

setwd("~/R-code/RNA-seq/No-Male-Enriched-Genes/")
source("Scripts/make.transparent.R")
library(Biobase)

files <- list.files("Data/Star-ReadCounts",
                    pattern = "ReadsPerGene",
                    full.names = TRUE)
## Remove male samples from the datasets
files <- files[c(1:28, 33:40)]

samples <- gsub("../Data/Star-ReadCounts/", "", files)
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
pData$bg <- make.transparent(pData$col, alpha = 200)
pData$shp[pData$Presence.of.males == "none"] <- 21
pData$shp[pData$Presence.of.males == "N2 males"
          & pData$Days.with.Males == "1"] <- 24
pData$shp[pData$Presence.of.males == "N2 males"
          & pData$Days.with.Males == "5"] <- 23
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
## 31041 15707
genes.keep <- names(ind.keep)[ind.keep]
eset <- eset[genes.keep, ]

table(fData(eset)$biotype)
##        antisense        lincRNA          miRNA          ncRNA          piRNA protein_coding     pseudogene 
##            12             63              6            128             18          15151            295
##        rRNA         snoRNA          snRNA           tRNA 
##        5             26              2               1 
dim(eset)
## Features  Samples
##   15707       36

## Coverage per sample:
quantile(colSums(exprs(eset)))
##       0%      25%      50%      75%     100%
##   5286965  9922946 12359191 18264467 25478068

## Number of genes per sample:
quantile(colSums(exprs(eset) != 0))
##      0%     25%     50%     75%    100%
##  14067.0 14487.5 14918.0 15238.0 15658.0

## Filter out genes that are expressed in more in males (N2 males) than glp-1 hermaphrodites (day 3 and 7 combined)
## The filtering is important because many of the most differentially expressed genes (+ vs. - male conditions) are sperm related
## In addition, some of the sample 'spread' or 'outliers' could be due to either a few males in the her samples, a non-glp-1 herm (technical problems)
## Another possibility is that some samples could have more mating (or more recent) and therefore more sperm
## The latter possibility is more interesting biologically and could provide hints for the role of mating and sperm/seminal fluid

load("Data/male-enriched-p05.RData")
male.genes <- sig.results.positive
eset.nomalegenes <- eset[!is.element(rownames(eset), male.genes), ]


table(fData(eset.nomalegenes)$biotype)
##        antisense        lincRNA          miRNA          ncRNA          piRNA protein_coding     pseudogene 
##             11             49              6            120             18           9923            19 
##        rRNA         snoRNA          snRNA           tRNA 
##          4             25              2              1
dim(eset.nomalegenes)
## Features  Samples
##  10352       36

## Coverage per sample:
quantile(colSums(exprs(eset.nomalegenes)))
##       0%      25%      50%      75%     100%
##    5072249  9485531 11858561 17187475 24474331 

## Number of genes per sample:
quantile(colSums(exprs(eset.nomalegenes) != 0))
##      0%     25%     50%     75%    100%
##    9856.00 10009.00 10093.50 10172.25 10318.00

eset <- eset.nomalegenes
save(eset,
     file = "Data/eset_counts.RData")

