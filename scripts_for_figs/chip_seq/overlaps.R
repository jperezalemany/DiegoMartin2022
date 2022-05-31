#!/usr/bin/env Rscript

# Venn diagrams of intersection between ChIP replicate peaks

setwd("~/../Dropbox/Bioinformatics/ChIP-seq TPF1 TPF2 CHR23/")

library(eulerr)
library(RColorBrewer)

colors <- brewer.pal(8, "Dark2")

tpf1 <- read.table("data/TPF1_peaks.tab", sep = "\t", header = T)
svg(file = "figures/overlap_replicates/TPF1.svg", height = 2, width = 2)
plot(
  venn(
    combinations = list(
      TPF1_rep1 = tpf1[tpf1$TPF1.4 == "True", ]$name,
      TPF1_rep2 = tpf1[tpf1$TPF1.11 == "True", ]$name
    )
  ),
  edges = list(col = colors[2:3], lex = 1.5),
  quantities = list(type = "counts", fontsize = 11),
  fills = NA
)
dev.off()

tpf2 <- read.table("data/TPF2_peaks.tab", sep = "\t", header = T)

svg(file = "figures/overlap_replicates/TPF2.svg", height = 2, width = 2)
plot(
  venn(
    combinations = list(
      TPF2_rep1 = tpf2[tpf2$TPF2.9 == "True", ]$name,
      TPF2_rep2 = tpf2[tpf2$TPF2.15 == "True", ]$name
    )
  ),
  edges = list(col = colors[2:3], lex = 1.5),
  quantities = list(type = c("counts"), fontsize = 11),
  fills = NA
)
dev.off()

chr23 <- read.table("data/CHR23_2_peaks.tab", sep = "\t", header = T)
head(chr23)

svg(file = "figures/overlap_replicates/CHR23.svg", height = 2, width = 2)
plot(
  venn(
    combinations = list(
      CHR23_rep1 = chr23[chr23$CHR23.7_2 == "True", ]$name,
      CHR23_rep2 = chr23[chr23$CHR23.10_2 == "True", ]$name
    )
  ),
  edges = list(col = colors[2:3], lex = 1.5),
  quantities = list(type = c("counts"), fontsize = 11),
  fills = NA
)
dev.off()

# Overlap between unionset of each protein
overlap <- read.table("data/TPF1_TPF2_CHR23_peaks.tab", sep = "\t", header = T)

svg(file = "figures/overlap_replicates/overlap_peaks.svg", height = 2, width = 2)
plot(
  venn(
    combination = list(
      TPF1 = overlap[overlap$TPF1 == "True", ]$name,
      TPF2 = overlap[overlap$TPF2 == "True", ]$name,
      CHR23 = overlap[overlap$CHR23 == "True", ]$name
    )
  ),
  edges = list(col = c(colors[1], colors[3], colors[2]), lex = 1.5),
  quantities = list(type = c("counts"), fontfamily = c("arial"), fontsize = c(11)),
  fills = NA
)
dev.off()

