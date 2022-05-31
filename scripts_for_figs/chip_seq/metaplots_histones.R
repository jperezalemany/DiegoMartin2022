#!/usr/bin/env Rscript

# TPF targets over histone marks
# Jaime Perez Alemany
# 07/05/2022

setwd("~/../Dropbox/Bioinformatics/ChIP-seq TPF1 TPF2 CHR23/")

library(dplyr)
library(cowplot)

profiles <- read.table("data/epigenetic_profile_2000.inter.tab", skip = 2)
dim(profiles)
colnames(profiles) <- c("chip", "genes", c(1:(ncol(profiles)-2)))

chips <- c("TPF1", "TPF2", "MINU2", "H3K4me3", "H3K27me3")

profiles[1:6, 1:6]
data <- profiles %>%
  pivot_longer(!c(chip, genes), names_to = "bin", values_to = "value") %>%
  mutate(genes = if_else(genes == "TPF.inter_targets.bed", "Targeted by TPF", "Not targeted by TPF")) %>%
  mutate(bin = as.integer(bin)) %>%
  mutate(chip = str_remove(chip, ".log2ratio")) %>%
  mutate(chip = if_else(chip == "TPF1_over_Col0", "TPF1", chip)) %>%
  mutate(chip = if_else(chip == "TPF2_over_Col0_1", "TPF2", chip)) %>%
  mutate(chip = if_else(chip == "CHR23_2_over_Col0_2", "MINU2", chip)) %>%
  mutate(chip = if_else(chip == "Col0_H3K4me3_over_Col0_H3", "H3K4me3", chip)) %>%
  mutate(chip = if_else(chip == "Col0_H3K27me3_over_Col0_H3", "H3K27me3", chip)) %>%
  filter(chip %in% chips) %>%
  mutate(chip = factor(chip, levels = chips))

plot1 <- data %>%
  filter(chip == "H3K4me3") %>%
  ggplot(aes(bin, value)) +
  geom_line(aes(color = genes), size = 1) +
  ylab(expression(Log[2]~"H3K4me3 / H3")) +
  theme_bw() +
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(2, 1)]) +
  scale_y_continuous(limits = c(-1, 2)) +
  scale_x_continuous(
    breaks = c(0, 200, 400, 600), labels = c("-2 kb", "TSS", "TTS", "2 kb")) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))
plot1
plot2 <- data %>%
  filter(chip == "H3K27me3") %>%
  ggplot(aes(bin, value)) +
  geom_line(aes(color = genes), size = 1) +
  ylab(expression(Log[2]~"H3K27me3 / H3")) +
  theme_bw() +
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(2, 1)]) +
  scale_y_continuous(limits = c(-0.7, 0.3), breaks = c(-0.6, -0.3, 0, 0.3)) +
  scale_x_continuous(
    breaks = c(0, 200, 400, 600), labels = c("-2 kb", "TSS", "TTS", "2 kb")) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))
  
plot2
plot_grid(plot1, plot2, ncol = 1, align = "v", axis = "rl")

ggsave("figures/metaplots/histones_template.pdf", dpi = 1200,
       height = 8, width = 6, units = "cm")
