#!/usr/bin/env Rscript

# Expression boxplots TPF targets vs non targets
# Jaime Perez Alemany
# 15/05/2022

library(dplyr)
library(ggplot2)
library(DESeq2)

setwd("~/../Dropbox/Bioinformatics/ChIP-seq TPF1 TPF2 CHR23/")
tpf_targets <- read.table("bedfiles/TPF.inter_targets.bed")
tpf_notargets <- read.table("bedfiles/TPF.inter_notargets.bed")

#dds_salmon <- readRDS("../RNA-seq tpf chr/deseq2/dds_salmon.RDS")
dds_counts <- readRDS("../RNA-seq tpf chr/deseq2/dds.RDS")
mcols(dds_counts)$basepairs <- mcols(dds_counts)$basepair

mrna_genes <- read.table("../Araport11/genes_mrna.bed")

#exp_salmon <- fpkm(dds_salmon) %>% data.frame() %>%
#  mutate(mean_col0 = (col0_rep1 + col0_rep2 + col0_rep3) / 3) %>%
#  tibble::rownames_to_column(var = "TAIR") %>%
#  select(TAIR, mean_col0) %>%
#  mutate(target = if_else(TAIR %in% tpf_targets$V4, "TPF-occ.",
#                          if_else(TAIR %in% tpf_notargets$V4, "TPF-free", "NA"))) %>%
#  mutate(mrna = if_else(TAIR %in% mrna_genes$V4, T, F))

exp_counts <- fpkm(dds_counts) %>% data.frame() %>%
  mutate(mean_col0 = (col0_rep1 + col0_rep2 + col0_rep3) / 3) %>%
  tibble::rownames_to_column(var = "TAIR") %>%
  select(TAIR, mean_col0) %>%
  mutate(target = if_else(TAIR %in% tpf_targets$V4, "TPF-occ.",
                          if_else(TAIR %in% tpf_notargets$V4, "TPF-free", "NA"))) %>%
  mutate(mrna = if_else(TAIR %in% mrna_genes$V4, T, F))

data <- exp_counts %>%
  filter(target != "NA") %>%
  mutate(target = factor(target, levels = c("TPF-occ.", "TPF-free"))) %>%
  filter(mrna == T) %>%
  mutate(log2_mean = log2(mean_col0)) %>%
  filter(!is.na(log2_mean))

plot <- exp_counts %>%
  filter(target != "NA") %>%
  mutate(target = factor(target, levels = c("TPF-occ.", "TPF-free"))) %>%
  ggplot(aes(target, log2(mean_col0), color = target)) +
  ylab(expression(Log[2]~FPKM)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  
  geom_boxplot(
    outlier.shape = 1, outlier.size = 1, width = 0.5) +
  scale_y_continuous(limits = c(-10, 15)) +
  theme_bw() +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11)
  )

dir.create("figures/expression/")
ggsave(plot, filename = "figures/expression/boxplot.pdf",
       height = 4, width = 5.5, units = "cm", dpi = 1200)


