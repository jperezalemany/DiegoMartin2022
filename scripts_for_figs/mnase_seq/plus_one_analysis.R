#!/usr/bin/env Rscript

# Shift analysis
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ggtext)

nuc_columns <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "smt_score", "smt_pval", "fuzz_score", "fuzz_pval",
  "gene_name", "tss_distance"
)

setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

greylist <- read.table("bedfiles/col0.merged_annotated_greylist.bed")

targets <- read.table("../ChIP-seq TPF1 TPF2 CHR23/bedfiles/TPF.inter_targets.bed")
notargets <- read.table("../ChIP-seq TPF1 TPF2 CHR23/bedfiles/TPF.inter_notargets.bed")

# Read shift table
shift <- read.csv("processed_tables/chr_tpf.integrated_shift.csv")
shift$col0_comb <- paste(shift$control_name, shift$gene_name)
shift$tpf_comb <- paste(shift$tpf_name, shift$gene_name)
shift$chr_comb <- paste(shift$chr_name, shift$gene_name)


# Plus ones
plus_one_col0 <- read.table("bedfiles/col0_plus_one.bed")
nrow(plus_one_col0)
plus_one_col0$comb <- paste(plus_one_col0$V4, plus_one_col0$V11)
colnames(plus_one_col0) <- paste0("col0_", c(nuc_columns, "comb"))

plus_one_chr <- read.table("bedfiles/chr_plus_one.bed")
plus_one_chr$comb <- paste(plus_one_chr$V4, plus_one_chr$V11)
colnames(plus_one_chr) <- paste0("chr_", c(nuc_columns, "comb"))

plus_one_tpf <- read.table("bedfiles/tpf_plus_one.bed")
plus_one_tpf$comb <- paste(plus_one_tpf$V4, plus_one_tpf$V11)
colnames(plus_one_tpf) <- paste0("tpf_", c(nuc_columns, "comb"))

# Merge
tpf_summits <- read.table("bedfiles/tpf_summits.bed")
colnames(tpf_summits) <- paste0("tpf_", nuc_columns[1:10])

chr_summits <- read.table("bedfiles/chr_summits.bed")
colnames(chr_summits) <- paste0("chr_", nuc_columns[1:10])

# Merge shift with plus one positions in WT
plus_one_data <- merge(plus_one_col0, shift, by = "col0_comb") %>%
  merge(tpf_summits, by = "tpf_name", all.x = T) %>%
  merge(chr_summits, by = "chr_name", all.x = T) %>%
  mutate(target = if_else(gene_name %in% targets$V4, T,
                           if_else(gene_name %in% notargets$V4, F, NA))) %>%
  arrange((chr2control_dis + tpf2control_dis)/2)

write.csv(plus_one_data, file = "processed_tables/plus_one_data.csv")

# Nucleosomes W to W

targets %>% nrow()
dir.create("plus_one")
plus_one_data %>%
  filter(target == T) %>%
  select(col0_chrom, col0_start, col0_end, col0_name, col0_score, col0_strand) %>%
  write.table("plus_one/all_targets.bed", sep = "\t", quote = F, col.names = F, row.names = F)

plus_one_data %>%
  filter(target == F) %>%
  select(col0_chrom, col0_start, col0_end, col0_name, col0_score, col0_strand) %>%
  write.table("plus_one/all_notargets.bed", sep = "\t", quote = F, col.names = F, row.names = F)

plus_one_data %>%
  filter(col0_fuzz_pval < log10(0.05)) %>%
  filter(target == T) %>%
  select(col0_chrom, col0_start, col0_end, col0_name, col0_score, col0_strand) %>%
  write.table("plus_one/wtfuzz_targets.bed", sep = "\t", quote = F, col.names = F, row.names = F)

plus_one_data %>%
  filter(col0_fuzz_pval < log10(0.05)) %>%
  filter(target == F) %>%
  select(col0_chrom, col0_start, col0_end, col0_name, col0_score, col0_strand) %>%
  write.table("plus_one/wtfuzz_notargets.bed", sep = "\t", quote = F, col.names = F, row.names = F)

plus_one_data %>%
  filter(col0_fuzz_pval < log10(0.05)) %>%
  filter(
    chr_fuzz_pval < log10(0.05) & tpf_fuzz_pval < log10(0.05)
  ) %>%
  filter(
    chr_smt_pval < log10(0.05) & tpf_smt_pval < log10(0.05)
  ) %>%
  filter(target == T) %>%
  select(col0_chrom, col0_start, col0_end, col0_name, col0_score, col0_strand) %>%
  write.table("plus_one/wtmfuzz_targets.bed", sep = "\t", quote = F, col.names = F, row.names = F)

plus_one_data %>%
  filter(col0_fuzz_pval < log10(0.05)) %>%
  filter(
    chr_fuzz_pval < log10(0.05) & tpf_fuzz_pval < log10(0.05)
  ) %>%
  filter(
    chr_smt_pval < log10(0.05) & tpf_smt_pval < log10(0.05)
  ) %>%
  filter(target == F) %>%
  select(col0_chrom, col0_start, col0_end, col0_name, col0_score, col0_strand) %>%
  write.table("plus_one/wtmfuzz_notargets.bed", sep = "\t", quote = F, col.names = F, row.names = F)

# Shift deciles
plus_one_data %>%
  filter(col0_fuzz_pval < log10(0.05)) %>%
  filter(
    chr_fuzz_pval < log10(0.05) & tpf_fuzz_pval < log10(0.05)
  ) %>%
  filter(
    chr_smt_pval < log10(0.05) & tpf_smt_pval < log10(0.05)
  ) %>%
  filter(target == T) %>%
  nrow()


# Violin plots

chr_plot <- plus_one_data %>%
  filter(col0_fuzz_pval < log10(0.05)) %>%
  filter(
    chr_fuzz_pval < log10(0.05) & tpf_fuzz_pval < log10(0.05)
  ) %>%
  filter(
    chr_smt_pval < log10(0.05) & tpf_smt_pval < log10(0.05)
  ) %>%
  pivot_longer(cols = c(chr2control_dis, tpf2control_dis), names_to = "mutant", values_to = "shift") %>%
  filter(mutant == "chr2control_dis") %>%
  filter(!(is.na(target))) %>%
  mutate(target = if_else(target == T, "TPF-occ.", "TPF-free")) %>%
  ggplot(aes(target, shift, color = target)) +
  facet_wrap(~mutant) +
  ylab("+1 position (*minu* - Col-0)") +
  coord_flip() +
  scale_y_continuous(limits = c(-100, 100),
                     breaks = c(-100, -50, 0, 50, 100),
                     minor_breaks = seq(-100, 100, by = 10)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[2:1]) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_markdown(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))

chr_plot
tpf_plot <- plus_one_data %>%
  filter(col0_fuzz_pval < log10(0.05)) %>%
  filter(
    chr_fuzz_pval < log10(0.05) & tpf_fuzz_pval < log10(0.05)
  ) %>%
  filter(
    chr_smt_pval < log10(0.05) & tpf_smt_pval < log10(0.05)
  ) %>%
  pivot_longer(cols = c(chr2control_dis, tpf2control_dis), names_to = "mutant", values_to = "shift") %>%
  filter(mutant == "tpf2control_dis") %>%
  filter(!(is.na(target))) %>%
  mutate(target = if_else(target == T, "TPF-occ.", "TPF-free")) %>%
  ggplot(aes(target, shift, color = target)) +
  facet_wrap(~mutant) +
  ylab("+1 position (*tpf* - Col-0)") +
  coord_flip() +
  scale_y_continuous(limits = c(-100, 100),
                     breaks = c(-100, -50, 0, 50, 100),
                     minor_breaks = seq(-100, 100, by = 10)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[2:1]) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_markdown(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))

tpf_plot
dir.create("figures/violin_plots")
ggsave(tpf_plot, filename = "figures/violin_plots/tpf_violin.pdf", units = "cm",
       height = 4, width = 8, dpi = 1200)
ggsave(chr_plot, filename = "figures/violin_plots/chr_violin.pdf", units = "cm",
       height = 4, width = 8, dpi = 1200)


plus_one_data %>%
  filter(target == T) %>%
  filter(mean_s)
plus_one_data <- read.csv("processed_tables/plus_one_data.csv")
# Average line

data %>%
  arrange(col0_fuzz_pval) %>%
  head()

data <- plus_one_data %>%
  filter(col0_fuzz_pval < log10(0.05)) %>%
  filter(
    chr_fuzz_pval < log10(0.05) & tpf_fuzz_pval < log10(0.05)
  ) %>%
  filter(
    chr_smt_pval < log10(0.05) & tpf_smt_pval < log10(0.05)
  ) %>%
  mutate(mean_shift = (chr2control_dis + tpf2control_dis) / 2) %>%
  filter(target == T)

data %>% filter(col0_strand == "+") %>% head(n = 200) %>% arrange(col0_fuzz_pval + chr_fuzz_pval + tpf_fuzz_pval) %>%
  head(n = 20)

data$rank <- c(1:nrow(data))

rank5 <- data[data$mean_shift == -5, "rank"][1]
rank10 <- data[data$mean_shift == -10, "rank"][1]
rank20 <- data[data$mean_shift == -20, "rank"][1]

rank5
ggplot(data, aes(-rank, mean_shift)) +
  scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(limits = c(-100, 100),
                     breaks = c(-100, -50, 0, 50, 100),
                     minor_breaks = c(-75, -25, 25, 75)) +
  geom_vline(xintercept = c(-rank5, -rank20, -rank10)) +
  geom_area() +
  xlab("+1 position avg.(*tpf*, *minu*) - WT") +
  coord_flip() +
  theme_bw() +
  theme(
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.line.x = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12)
  )

ggsave(filename = "figures/heatmaps/shift_line_test.pdf",
       width = 40, height = 53.921, units = "mm", dpi = 1200)
