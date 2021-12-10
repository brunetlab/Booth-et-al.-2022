## This script will create a balloon plot to sumarize my MID published data

rm(list=ls())
library(ggplot2)
setwd("~/Dropbox/MID-lifespan-analysis/Summary-Figures/Screen+longevity/")

results <- read.delim("Screen_results_cox_with_longevity.txt", 
                             stringsAsFactors = FALSE,
                             row.names = 1)

## bulble plot of all results:
a <- ggplot(results, aes(x = gene, y = condition, size = limited.neglog10pval, fill = limited.HR))+
  geom_point(shape = 21) +
  ggtitle("Screen Results and Longevity Pathways") +
  labs(x = "Gene") +
  theme(panel.background=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(7, "PRGn")),
                       breaks=c(0.001, 0.05, 0.5, 1, 1.5, 2, 2.5)) +
  scale_x_discrete(limits = results$gene[order(results$x.order)]) +
  theme(axis.text.x=element_text(angle=90, hjust=1))
a

pdf("summary_cox.pdf",
    width = 10, height = 3, 
    onefile=FALSE,
    useDingbats = FALSE)
plot(a)
dev.off()

