rm(list=ls())
library(ggplot2)
library(RColorBrewer)

setwd("~/Lifespan-Analysis/")

cox.interaction <- read.delim("Interaction-Results.txt")

## create box plot for the interaction results:
box.plot.basic <- ggplot(cox.interaction, aes(x = Gene, y = HR)) +
  geom_boxplot(aes(x=reorder(Gene,HR), y=HR), outlier.shape = NA, lwd = 0.4) +
  geom_point(aes(shape = Label, color = "grey", alpha = 70), position = position_jitterdodge()) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks=c(0.25, 0.5, 1, 10, 20, 30)) +
  labs(title = "Lifespan Data: Cox Proportional Hazard Ratios, gene x MID", x = "Gene", y = "Hazard Ratio")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("Figure2-boxplot-interaction.pdf",
    width = 9, height = 3, 
    onefile=FALSE,
    useDingbats = FALSE)
plot(box.plot.basic)

dev.off()


## create box plot for the +males and her only results:
cox.her.males <- read.delim("Single-Condition-Results.txt")

box.plot.basic <- ggplot(cox.her.males, aes(x = Gene, y = HR)) +
  geom_boxplot(aes(x=reorder(Gene,Ordering), y=HR, fill = Condition), outlier.shape = NA, lwd = 0.4) +
  geom_point(aes(shape = Label, color = Condition, alpha = 70), position = position_jitterdodge()) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  scale_y_continuous(breaks=c(0.01, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75)) +
  labs(title = "Lifespan Data: Cox Proportional Hazard Ratios", x = "Gene", y = "Hazard Ratio")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("Figure2-boxplot-single-condition.pdf",
    width = 10, height = 3, 
    onefile=FALSE,
    useDingbats = FALSE)
plot(box.plot.basic)

dev.off()

box.plot.basic <- ggplot(cox.her.males, aes(x = Gene, y = HR)) +
  geom_boxplot(aes(x=reorder(Gene,Ordering), y=HR, fill = Condition), outlier.shape = NA, lwd = 0.4) +
  geom_point(aes(shape = Label, color = Condition, alpha = 70), position = position_jitterdodge()) +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks=c(0.25, 0.5, 1, 10, 20, 30)) +
  labs(title = "Lifespan Data: Cox Proportional Hazard Ratios", x = "Gene", y = "Hazard Ratio")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("Figure2-boxplot-single-condition-log10.pdf",
    width = 10, height = 3, 
    onefile=FALSE,
    useDingbats = FALSE)
plot(box.plot.basic)

dev.off()
