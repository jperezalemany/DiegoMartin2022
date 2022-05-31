
library(DESeq2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("~/../Dropbox/Bioinformatics/RNA-seq tpf chr")
dds <- readRDS("deseq2/dds.RDS")

rld <- rlog(dds, blind = F)

plotPCA(rld, ntop = Inf)
pcadata <- plotPCA(rld, ntop = Inf, returnData = T)
pcadata$replicate <- rep(c(1, 2, 3), 3)
ggplot(pcadata, aes(PC1, PC2, color = condition)) +
  geom_point(size = 2) +
  xlab("PC1: 65% variance") +
  ylab("PC2: 19 % variance") +
  theme_bw() +
  scale_color_manual(values = brewer.pal(3, "Set2")) +
  scale_x_continuous(limits = c(-60, 60)) +
  scale_y_continuous(limits = c(-40, 40)) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11)
        ) +
  coord_fixed()

dir.create("figures/supp")
ggsave(filename = "figures/supp/pca.pdf", height = 6, width = 6, units = "cm", dpi = 1200)
