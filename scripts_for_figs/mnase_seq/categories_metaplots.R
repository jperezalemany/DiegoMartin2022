#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)

setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

data <- read.table("data/occ_500.tab", skip = 2, sep = "\t")
colnames(data) <- c("background", "genes", c(1:1000))

data <- data %>%
  pivot_longer(cols = as.character(c(1:1000)), names_to = "bin", values_to = "value") %>%
  mutate(bin = as.integer(bin))

tpf_occ <- data %>%
  filter(bin >= 250, bin <= 750) %>%
  filter(background == "col0" | background == "tpf") %>%
  filter(genes == "tpf_up" | genes == "tpf_down") %>%
  mutate(genes = factor(genes, levels = c("tpf_up", "tpf_down"))) %>%
  ggplot(aes(bin, value)) +
  geom_line(aes(color = background)) +
  scale_color_manual(values = c("grey", brewer.pal(3, "Dark2")[3])) +
  scale_x_continuous(breaks = c(250, 500, 750), labels = c(-250, 0, 250)) +
  scale_y_continuous(limits = c(2500, 7000), breaks = c(2500, 7000)) +
  ylab("Norm. read count") +
  xlab("Distance from dyad summit (Kb)") +
  facet_wrap(~genes) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(color = "black", size = 11))

chr_occ <- data %>%
  filter(bin >= 250, bin <= 750) %>%
  filter(background == "col0" | background == "chr") %>%
  filter(genes == "chr_up" | genes == "chr_down") %>%
  mutate(genes = factor(genes, levels = c("chr_up", "chr_down"))) %>%
  mutate(background = factor(background, levels = c("col0", "chr"))) %>%
  ggplot(aes(bin, value)) +
  geom_line(aes(color = background)) +
  scale_color_manual(values = c("grey", brewer.pal(3, "Dark2")[2])) +
  scale_x_continuous(breaks = c(250, 500, 750), labels = c(-250, 0, 250)) +
  scale_y_continuous(limits = c(2500, 7000), breaks = c(2500, 7000)) +
  ylab("Norm. read count") +
  xlab("Distance from dyad summit (Kb)") +
  facet_wrap(~genes) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(color = "black", size = 11))
chr_occ
ggsave(chr_occ, filename = "figures/categories/occ_chr_metaplot.pdf", height = 3, width = 6, units = "cm", dpi = 1200)
