#!/usr/bin/env Rscript

# Occupancy analysis

library(dplyr)
library(RColorBrewer)
library(ggplot2)

setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

occ <- read.csv("processed_tables/chr_tpf.integrated_occupancy.csv") %>%
  filter(!is.na(chr_smt_diff_FDR)) %>%
  filter(!is.na(tpf_smt_diff_FDR))

occ %>% write.table("processed_tables/chr_tpf.integrated_occupancy.tab", sep = "\t", quote = F, row.names = F)
col0_summits <- read.table("bedfiles/col0_summits.annotation.bed") %>%
  filter(V8 < log10(0.05))

tpf_summits <- read.table("bedfiles/tpf_summits.annotation.bed") %>%
  filter(V8 < log10(0.05))

chr_summits <- read.table("bedfiles/chr_summits.annotation.bed") %>%
  filter(V8 < log10(0.05))

tpf_nucs <- c(col0_summits$V4, occ[occ$tpf_name %in% tpf_summits$V4, "control_name"]) %>% unique()
chr_nucs <- c(col0_summits$V4, occ[occ$chr_name %in% chr_summits$V4, "control_name"]) %>% unique()

occ$cat_chr <- ""
occ[occ$chr_smt_diff_FDR < 0.05 & occ$chr_smt_log2FC < 0 & occ$control_name %in% chr_nucs, "cat_chr"] <- "down"
occ[occ$chr_smt_diff_FDR < 0.05 & occ$chr_smt_log2FC > 0 & occ$control_name %in% chr_nucs, "cat_chr"] <- "up"

occ$cat_tpf <- ""
occ[occ$tpf_smt_diff_FDR < 0.05 & occ$tpf_smt_log2FC < 0 & occ$control_name %in% tpf_nucs, "cat_tpf"] <- "down"
occ[occ$tpf_smt_diff_FDR < 0.05 & occ$tpf_smt_log2FC > 0 & occ$control_name %in% tpf_nucs, "cat_tpf"] <- "up"


occ_data <- occ %>%
  filter(control_name %in% c(tpf_nucs, chr_nucs))

cor(occ_data$chr_smt_log2FC, occ_data$tpf_smt_log2FC)
pdf("figures/categories/occ_fc.pdf", height = 4, width = 4)
smoothScatter(occ_data$tpf_smt_log2FC, occ_data$chr_smt_log2FC, xlim = c(-12, 12), ylim = c(-12, 12), nbin = 100)
dev.off()

max(abs(occ_data$tpf_smt_log2FC))

ggplot(occ_data, aes(tpf_smt_log2FC, chr_smt_log2FC), ) +
  geom_hex(bins = 100, aes(colour = ..count..), size = 0.05) +
  theme_bw() +
  xlab("Log<sub>2</sub> (*tpf* / Col-0)") +
  ylab("Log<sub>2</sub> (*minu* / Col-0)") +
  scale_x_continuous(limits = c(-15, 15), breaks = c(-15, -7.5, 0, 7.5, 15)) +
  scale_y_continuous(limits = c(-15, 15), breaks = c(-15, -7.5, 0, 7.5, 15)) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(color = "black", size = 11)) +
  coord_fixed() +
  scale_fill_viridis_c() +
  scale_color_viridis_c()

dir.create("figures/categories/")
ggsave(filename = "figures/categories/occ_fc.pdf", height = 6, width = 6, units = "cm", dpi = 1200)
ggsave(filename = "figures/categories/occ_fc.tiff", height = 6, width = 6, units = "cm", dpi = 1200)

# Write table
data.frame(
  mutant = c("tpf", "chr", "tpf", "chr"),
  group = c("increase", "increase", "decrease", "decrease"),
  number = c(
    length(unique(occ[occ$cat_tpf == "down", "control_name"])),
    length(unique(occ[occ$cat_chr == "down", "control_name"])),
    length(unique(occ[occ$cat_tpf == "up", "control_name"])),
    length(unique(occ[occ$cat_chr == "up", "control_name"]))
  ),
  from = c(length(tpf_nucs), length(chr_nucs), length(tpf_nucs), length(chr_nucs))
) %>%
  write.table("processed_tables/occ.txt", sep = "\t", row.names = F, quote = F)

# Export
col0_summits <- read.table("bedfiles/col0_summits.annotation.bed")

col0_summits %>%
  filter(V4 %in% unique(occ[occ$cat_chr == "down", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V5) %>%
  distinct() %>%
  write.table("categories/occ_chr_down.bed", sep = "\t",
              col.names = F, row.names = F, quote = F)

col0_summits %>%
  filter(V4 %in% unique(occ[occ$cat_chr == "up", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V5) %>%
  distinct() %>%
  write.table("categories/occ_chr_up.bed", sep = "\t",
              col.names = F, row.names = F, quote = F)

col0_summits %>%
  filter(V4 %in% unique(occ[occ$cat_tpf == "down", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V5) %>%
  distinct() %>%
  write.table("categories/occ_tpf_down.bed", sep = "\t",
              col.names = F, row.names = F, quote = F)

col0_summits %>%
  filter(V4 %in% unique(occ[occ$cat_tpf == "up", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V5) %>%
  distinct() %>%
  write.table("categories/occ_tpf_up.bed", sep = "\t",
              col.names = F, row.names = F, quote = F)


