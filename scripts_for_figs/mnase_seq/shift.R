#!/usr/bin/env Rscript

# Shift analysis
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(scales)

setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

# Read shift table
shift <- read.csv("processed_tables/chr_tpf.integrated_shift.csv") %>%
  filter(!is.na(chr2control_dis)) %>%
  filter(!is.na(tpf2control_dis))

shift %>%
  write.table("processed_tables/chr_tpf.integrated_shift.tab", sep = "\t", row.names = F, quote = F)

col0_summits <- read.table("bedfiles/col0_summits.annotation.bed") %>%
  filter(V10 < log10(0.05)) %>%
  filter(V8 < log10(0.05))

tpf_summits <- read.table("bedfiles/tpf_summits.annotation.bed") %>%
  filter(V10 < log10(0.05)) %>%
  filter(V8 < log10(0.05))

chr_summits <- read.table("bedfiles/chr_summits.annotation.bed") %>%
  filter(V10 < log10(0.05)) %>%
  filter(V8 < log10(0.05))

chr_nucs <- intersect(col0_summits$V4, shift[shift$chr_name %in% chr_summits$V4, "control_name"]) %>% unique()
tpf_nucs <- intersect(col0_summits$V4, shift[shift$tpf_name %in% tpf_summits$V4, "control_name"]) %>% unique()

shift$cat_chr <- ""
shift[shift$chr2control_dis < -5 & shift$control_name %in% chr_nucs, "cat_chr"] <- "down"
shift[shift$chr2control_dis > 5 & shift$control_name %in% chr_nucs, "cat_chr"] <- "up"

table(shift$cat_chr)

shift$cat_tpf <- ""
shift[shift$tpf2control_dis < -5 & shift$control_name %in% tpf_nucs, "cat_tpf"] <- "down"
shift[shift$tpf2control_dis > 5 & shift$control_name %in% tpf_nucs, "cat_tpf"] <- "up"

#
shift_data <- shift %>%
  filter(control_name %in% c(tpf_nucs, chr_nucs))

cor(shift_data$chr2control_dis, shift_data$tpf2control_dis)

pdf("figures/categories/shift.pdf", height = 4, width = 4)
smoothScatter(x = shift_data$tpf2control_dis, y = shift_data$chr2control_dis, nbin = 100)
dev.off()

ggplot(shift_data, aes(tpf2control_dis, chr2control_dis)) +
  geom_hex(bins = 100, aes(colour = ..count..), size = 0.05) +
  theme_bw() +
  xlab("*tpf* - Col-0") +
  ylab("*minu* - Col-0") +
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
ggsave(filename = "figures/categories/shift.pdf", height = 6, width = 6, units = "cm", dpi = 1200)
ggsave(filename = "figures/categories/shift.tiff", height = 6, width = 6, units = "cm", dpi = 1200)


fuzz_data %>% filter(cat_chr == "up" & cat_tpf == "down") %>% arrange(tpf_fuzziness_diff_FDR + chr_fuzziness_diff_FDR) %>% head()

occ_data %>% filter(cat_chr == "down" & cat_tpf == "down") %>% arrange(chr_smt_diff_FDR + tpf_smt_diff_FDR) %>% head()


shift_data %>% filter(cat_chr == "down" & cat_tpf == "down") %>% 
  arrange(chr2control_dis + tpf2control_dis) %>% 
  head(n = 100)
# Write table

data.frame(
  mutant = c("tpf", "chr", "tpf", "chr"),
  group = c("upstream", "upstream", "downstream", "downstream"),
  number = c(
    length(unique(shift[shift$cat_tpf == "down", "control_name"])),
    length(unique(shift[shift$cat_chr == "down", "control_name"])),
    length(unique(shift[shift$cat_tpf == "up", "control_name"])),
    length(unique(shift[shift$cat_chr == "up", "control_name"]))
  ),
  from = c(length(tpf_nucs), length(chr_nucs), length(tpf_nucs), length(chr_nucs))
) %>%
  write.table("processed_tables/shift.txt", sep = "\t", row.names = F, quote = F)


# Export
col0_summits <- read.table("bedfiles/col0_summits.annotation.bed")

col0_summits %>%
  filter(V4 %in% unique(shift[shift$cat_chr == "down", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V6) %>%
  write.table("categories/shift_chr_down.bed", quote = F, col.names = F,
              row.names = F, sep = "\t")

col0_summits %>%
  filter(V4 %in% unique(shift[shift$cat_chr == "up", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V6) %>%
  write.table("categories/shift_chr_up.bed", quote = F, col.names = F,
              row.names = F, sep = "\t")

col0_summits %>%
  filter(V4 %in% unique(shift[shift$cat_tpf == "down", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V6) %>%
  write.table("categories/shift_tpf_down.bed", quote = F, col.names = F,
              row.names = F, sep = "\t")

col0_summits %>%
  filter(V4 %in% unique(shift[shift$cat_tpf == "up", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V6) %>%
  write.table("categories/shift_tpf_up.bed", quote = F, col.names = F,
              row.names = F, sep = "\t")
