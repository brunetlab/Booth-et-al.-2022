## This script will create balloon plots to sumarize my MID screen data

rm(list=ls())
library(ggplot2)
setwd("~/R-code/Lifespan-BalloonPlot/")

screen.results <- read.delim("screen-results.txt", 
                             stringsAsFactors = FALSE,
                             row.names = 1)


## Creates a balloon plot using the difference in maximum lifespan between control and RNAi treatment
## Presented in Figure 1i
p <- ggplot(screen.results, aes(x = gene, y = condition, size = neglog10pval, fill = max.age)) +
  geom_point(shape = 21) +
  ggtitle("Screen Results") +
  labs(x = "Gene") +
  theme(panel.background=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_gradient2(low = "purple4", high = "darkorange1", mid = "white", midpoint = 0) +
  scale_x_discrete(limits = screen.results$gene[order(screen.results$max.order)]) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
p

pdf("screen_summary_max.pdf",
    width = 10, height = 3, 
    onefile=FALSE,
    useDingbats = FALSE)
plot(p)

dev.off()


## Creates a balloon plot using the difference in median lifespan between control and RNAi treatment
## Presented in Supplemental Figure 2a
b <- ggplot(screen.results, aes(x = gene, y = condition, size = neglog10pval, fill = median.age)) +
  geom_point(shape = 21) +
  ggtitle("Screen Results") +
  labs(x = "Gene") +
  theme(panel.background=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_gradient2(low = "purple4", high = "darkorange1", mid = "white", midpoint = 0) +
  scale_x_discrete(limits = screen.results$gene[order(-screen.results$median.order)]) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
b

pdf("screen_summary_median.pdf",
    width = 10, height = 3, 
    onefile=FALSE,
    useDingbats = FALSE)
plot(b)

dev.off()