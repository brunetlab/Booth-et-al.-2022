Here I am only analyzing the samples that were in a hermaphrodite only setting (all ages) and the hermaphrodites exposed to N2 males (all ages). The transcripts that I defined as N2 male enriched (see next folder up) are excluded here to focus on genes changing in the hermaphrodite in response to males.  


1. Import.R
This script takes the read count values, sample information, and gene information to organize the data for DEseq2 and calculates basic information about the datasets (coverage, gene type counts etc.).  Import.R filters the data so that genes with low coverage are removed as well as the male-enriched genes. To focus on the genes that are changing in the hermaphrodites and not genes that are expressed in male sperm and detected in the hermaphrodites, I removed the male-enriched genes (adjusted p-value <= 0.05 and fold change > 0). These male-enriched genes were identified int eh next folder up. This filtered dataset was carried through the normal RNA-seq pipeline and was used for all downstream analysis here. 

2. Normalization.R
This script uses DEseq2 to calculate differential expression and generate a data dispersion plot for data quality testing. Normalize.R normalizes the datasets using variance stabilization transformation. This normalized data is used to generate PCA plots. Three sets of normalization were performed: all the samples in this pipeline, only the day3 hermaphrodite samples, and only the day 7 hermaphrodite samples

3. Descriptive.R and Descriptive-day*.R
These scripts perform principal component analysis and generates PCA plots and a dendrogram of the data. Though only one plot per script is presented in the paper, there are several PCA plots created with these scripts that show different PCs and color code the samples in different ways (for example, by library batch). Note that the colors are different and were modified in illustrator. 

4. Model_with_Groups.R
These scripts use DEseq2 to calculate differential expression. The values calculated here are used in GO_term_enrichment_*.R and are included in a supplemental table. 

5. Strip_Plot.R
This creates a jitter plot of the log2(fold change) for the +males vs. no males results and color codes the points based on their adjusted p-values.

6. Heatmap-log2FC-intersection.R
This script generates a heat map using log2(fold change) and adjusted p-values from DEseq2. These values are also used to determine which genes to use to build the heat map (based on whether they are significantly changing in more than one comparison).

7. Intersection-w-self-spermRNA-seq.R
Here, I compare the lists of genes that are differentially expressed in this RNA-seq experiment (using glp-1 sterile hermaphrodites and WT N2 males) and in our Self-sperm paper (Booth et al. eLife 2019) RNA-seq experiments (using WT N2 and fem-1 feminized worms and him-5 males). This script creates heat maps and gene lists. The heat maps are of the log2(fold change) numbers

8. Pearson-MF-SS-RNA-seq.R
Here, I am comparing the different RNA-seq data sets from this paper and the Booth et al. 2019 eLife paper. This script calculates Pearson's correlation between the DEseq2 results (+males vs. no males) and creates a heat map of the results.



