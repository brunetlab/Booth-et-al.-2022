## This builds a strip plot with the significant genes colored

setwd("~/R-code/RNA-seq/No-Male-Enriched-Genes/")
library(Biobase)

load("Data/males_v_nomales_deseq.RData")
age.results.N2 <- list(N2malesvctrl.3 = N2malesvctrl.3, 
                    N2malesvctrl.7.1 = N2malesvctrl.7.1, 
                    N2malesvctrl.7.5 = N2malesvctrl.7.5)

## Remove rows with NA p-values:
age.results.N2 <- lapply(age.results.N2,
                      function(x) {
                        x[!is.na(x$padj), ]
                      })
n <- sapply(age.results.N2, nrow)
names(n) <- names(age.results.N2)

## Order by pvalue:
age.results.N2 <- lapply(age.results.N2,
                      function(x) {
                        x[order(x$padj),]
                      })

p.sig <- 0.05
cols <- list()
xlab <- character(length = length(age.results.N2))
for(i in seq(along = age.results.N2)){
  cols[[i]] <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 40), n[i]) # grey60
  ind.sig.i <- age.results.N2[[i]]$padj < p.sig
  cols[[i]][ind.sig.i] <- rgb(255, 0, 0, alpha=40, maxColorValue = 255)
  xlab[i] <- paste(names(age.results.N2)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(age.results.N2)
legend <- paste("p <", "", p.sig)
legend <- c(paste("p >", "", p.sig), legend)

pdf("Results/stripplot_deseq_Y-v-O-her+N2males-p<0.05.pdf", width = 5, height = 5)
par(mar = c(8.1, 4.1, 1, 7.1),
    xpd = TRUE)
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 3.5),
     ylim = c(-8, 8),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change (+WT males / no males)"
)
abline(h = 0,
       xpd = FALSE)
abline(h = seq(-8, 8, by = 2),
       lty = "dotted",
       col = "grey",
       xpd = FALSE)
for(i in 1:length(age.results.N2)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results.N2[[i]])), amount = 0.3),
         y = rev(age.results.N2[[i]]$log2FoldChange),
         pch = 21,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]))
}
axis(1,
     at = 1:3,
     tick = FALSE,
     las = 3,
     lwd = 0,
     labels = xlab)
axis(2,
     las = 1,
     at = seq(-8, 8, 2))
legend(x = 3.8, y = 5,
       pch = 21,
       col = unique(rev(cols[[i]])),
       pt.bg = unique(rev(cols[[i]])),
       bty = "n",
       legend = legend)
dev.off()

