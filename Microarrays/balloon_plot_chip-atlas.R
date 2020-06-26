## This script will create a balloon plot to summarize selected ChIP Atlas results

rm(list=ls())
library(ggplot2)
setwd("~/R-code/Microarrays/")

chip.results <- read.delim("Data/Results-2kb-sig.txt", 
                             stringsAsFactors = FALSE,
                             row.names = 1)

p <- ggplot(chip.results, aes(x = DEG, y = Antigen, size = enrichment.maxed, fill = Significance)) +
  geom_point(shape = 21) +
  ggtitle("Transcription factor binding enrichment") +
  labs(y = "Transcription factor") +
  ##scale_size_binned() +
  theme_bw()+
  ##theme(panel.background=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_gradient2(low = "deeppink3", mid = "yellow", midpoint = -1, na.value = "grey") +
  scale_y_discrete(limits = chip.results$Antigen[order(-chip.results$y.order)], expand = c(.05, 0)) +
  scale_x_discrete(limits = chip.results$DEG[order(chip.results$x.order)], expand = c(.1, 0)) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
p

pdf("summary.pdf",
    width = 5, height = 6, 
    onefile=FALSE,
    useDingbats = FALSE)
plot(p)

dev.off()