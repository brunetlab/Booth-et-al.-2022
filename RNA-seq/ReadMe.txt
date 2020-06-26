Here I am analyzing the samples that were in a hermaphrodite only setting (all ages), the N2 males, and the hermaphrodites exposed to N2 males (all ages). This data is also used to define a set of N2 male-enriched transcripts.

The scripts for this set of data analysis are:

1. Import.R
This script takes the read count values, sample information, and gene information to organize the data for DEseq2. Import.R filters the data so that genes with low coverage are removed. This script also calculates basic information about the datasets (coverage, gene type counts etc.). 

2. Normalization.R
This script uses DEseq2 to calculate differential expression and generate a data dispersion plot for data quality testing. Normalize.R normalizes the datasets using variance stabilization transformation. This normalized data is used to generate PCA plots. 

3. Descriptive.R
This script performs principal component analysis and generates PCA plots and a dendrogram of the data. There are several PCA plots created that show different PCs and color code the samples in different ways (for example, by library batch).

4. Model_with_Groups.R and Model_with_Male-v-Her_Groups.R
These scripts use DEseq2 to calculate differential expression. The values calculated here are used in GO_term_enrichment_*.R and are included in a supplemental table. Two types of groups are used: Model_with_Groups.R compares samples on the basis of sex, age, and exposure to the opposite sex whereas Model_with_Male-v-Her_Groups.R only compares males to hermaphrodites that never interacted with males. The latter is used to define a set of male-enriched genes that are used to filter the data in the next set of analyses.

(5. make.transparent.R is only used for aesthetics)


