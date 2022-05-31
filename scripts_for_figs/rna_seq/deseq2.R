#!/usr/bin/env Rscript

# tpf1 tpf2 and minu1 minu2 DEGs
# 12/05/2022

setwd("~/../Dropbox/Bioinformatics/RNA-seq tpf chr/")

library(DESeq2)
library(ggplot2)
library(dplyr)

dds <- readRDS("deseq2/dds.RDS")

# PCA figure
rld <- rlog(dds, blind = F)

plotPCA(rld, intgroup = "condition", ntop = Inf)
ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point(size = 1)

# Results
tpf <- lfcShrink(dds, coef = "condition_tpf_vs_col0")
chr <- lfcShrink(dds, coef = "condition_chr_vs_col0")

# Data frames
tpf_df <- tpf %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "TAIR") %>%
  mutate(deg = if_else(padj < 0.05 & log2FoldChange < -log2(1.5), "Down",
                       if_else(padj < 0.05 & log2FoldChange > log2(1.5), "Up", "NS")
                               )
                       )
table(tpf_df$deg)

chr_df <- chr %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "TAIR") %>%
  mutate(deg = if_else(padj < 0.05 & log2FoldChange < -log2(1.5), "Down",
                       if_else(padj < 0.05 & log2FoldChange > log2(1.5), "Up", "NS"))
  )

table(chr_df$deg)
# Save
write.csv(tpf_df, file = "deseq2/tpf_vs_col0.counts.csv", row.names = F, quote = F)
write.csv(chr_df, file = "deseq2/chr_vs_col0.counts.csv", row.names = F, quote = F)

