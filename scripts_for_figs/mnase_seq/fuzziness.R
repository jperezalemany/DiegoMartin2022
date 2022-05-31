#!/usr/bin/env Rscript

setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

library(dplyr)
library(ggtext)
library(ggplot2)

fuzz <- read.csv("processed_tables/chr_tpf.integrated_fuzziness.csv") %>%
  filter(!is.na(chr_fuzziness_diff_FDR)) %>%
  filter(!is.na(tpf_fuzziness_diff_FDR))

fuzz %>%
  write.table("processed_tables/chr_tpf.integrated_fuzziness.tab", sep = "\t", quote = F, row.names = F)

col0_summits <- read.table("bedfiles/col0_summits.annotation.bed") %>%
  filter(V8 < log10(0.05))

tpf_summits <- read.table("bedfiles/tpf_summits.annotation.bed") %>%
  filter(V8 < log10(0.05))

chr_summits <- read.table("bedfiles/chr_summits.annotation.bed") %>%
  filter(V8 < log10(0.05))

tpf_nucs <- c(col0_summits$V4, fuzz[fuzz$tpf_name %in% tpf_summits$V4, "control_name"]) %>% unique()
chr_nucs <- c(col0_summits$V4, fuzz[fuzz$chr_name %in% chr_summits$V4, "control_name"]) %>% unique()

fuzz$cat_chr <- ""
fuzz[fuzz$chr_fuzziness_diff_FDR < 0.05 & fuzz$chr_fuzziness_log2FC < 0 & fuzz$control_name %in% chr_nucs, "cat_chr"] <- "down"
fuzz[fuzz$chr_fuzziness_diff_FDR < 0.05 & fuzz$chr_fuzziness_log2FC > 0 & fuzz$control_name %in% chr_nucs, "cat_chr"] <- "up"

fuzz$cat_tpf <- ""
fuzz[fuzz$tpf_fuzziness_diff_FDR < 0.05 & fuzz$tpf_fuzziness_log2FC < 0 & fuzz$control_name %in% tpf_nucs, "cat_tpf"] <- "down"
fuzz[fuzz$tpf_fuzziness_diff_FDR < 0.05 & fuzz$tpf_fuzziness_log2FC > 0 & fuzz$control_name %in% tpf_nucs, "cat_tpf"] <- "up"

# Fig
table(fuzz$cat_chr)

fuzz_data <- fuzz %>%
  filter(control_name %in% c(tpf_nucs, chr_nucs))
cor(fuzz_data$chr_fuzziness_log2FC, fuzz_data$tpf_fuzziness_log2FC)

pdf("figures/categories/fuzz_fc.pdf", height = 4, width = 4)
smoothScatter(x = fuzz_data$tpf_fuzziness_log2FC, y = fuzz_data$chr_fuzziness_log2FC, xlim = c(-1, 1), ylim = c(-1, 1), nbin = 100)
dev.off()

smoo

max(abs(fuzz_data$tpf_fuzziness_log2FC))
ggplot(fuzz_data, aes(tpf_fuzziness_log2FC, chr_fuzziness_log2FC)) +
  xlim(-1.1, 1.1) +
  ylim(-1.1, 1.1) +
  geom_hex(bins = 100, aes(colour = ..count..), size = 0.05) +
  
  theme_bw() +
  xlab("Log<sub>2</sub> (*tpf* / Col-0)") +
  ylab("Log<sub>2</sub> (*minu* / Col-0)") +
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
ggsave(filename = "figures/categories/fuzz_fc.pdf", height = 6, width = 6, units = "cm", dpi = 1200)
ggsave(filename = "figures/categories/fuzz_fc.tiff", height = 6, width = 6, units = "cm", dpi = 1200)


data.frame(
  mutant = c("tpf", "chr", "tpf", "chr"),
  group = c("increase", "increase", "decrease", "decrease"),
  number = c(
    length(unique(fuzz[fuzz$cat_tpf == "down", "control_name"])),
    length(unique(fuzz[fuzz$cat_chr == "down",  "control_name"])),
    length(unique(fuzz[fuzz$cat_tpf == "up",  "control_name"])),
    length(unique(fuzz[fuzz$cat_chr == "up",  "control_name"]))
  ),
  from = c(length(tpf_nucs), length(chr_nucs), length(tpf_nucs), length(chr_nucs))
) %>%
  write.table("processed_tables/fuzz.txt", sep = "\t", row.names = F, quote = F)

# Export
col0_summits <- read.table("bedfiles/col0_summits.annotation.bed")

col0_summits %>%
  filter(V4 %in% unique(fuzz[fuzz$cat_chr == "down", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V5) %>%
  distinct() %>%
  write.table("categories/fuzz_chr_down.bed", sep = "\t",
              col.names = F, row.names = F, quote = F)

col0_summits %>%
  filter(V4 %in% unique(fuzz[fuzz$cat_chr == "up", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V5) %>%
  distinct() %>%
  write.table("categories/fuzz_chr_up.bed", sep = "\t",
              col.names = F, row.names = F, quote = F)

col0_summits %>%
  filter(V4 %in% unique(fuzz[fuzz$cat_tpf == "down", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V5) %>%
  distinct() %>%
  write.table("categories/fuzz_tpf_down.bed", sep = "\t",
              col.names = F, row.names = F, quote = F)

col0_summits %>%
  filter(V4 %in% unique(fuzz[fuzz$cat_tpf == "up", "control_name"])) %>%
  select(V1, V2, V3, V4, V5, V5) %>%
  distinct() %>%
  write.table("categories/fuzz_tpf_up.bed", sep = "\t",
              col.names = F, row.names = F, quote = F)

    
