#!/usr/bin/Rscript

# TPF1 TPF2 CHR23 signal over deciles of expression
# Jaime Perez Alemany
# 15/05/2022

library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(RColorBrewer)

setwd("~/../Dropbox/Bioinformatics/ChIP-seq TPF1 TPF2 CHR23/")

profiles <- read.table("data/fpkm_deciles_counts.tab", skip = 2)
colnames(profiles) <- c("chip", "decile", c(1:600))

data <- profiles %>%
  pivot_longer(as.character(c(1:600)), names_to = "bin", values_to = "signal") %>%
  mutate(bin = as.integer(bin)) %>%
  mutate(decile = as.factor(decile)) %>%
  mutate(chip = factor(chip, levels = c("TPF1", "TPF2", "CHR23")))

fig <- data %>%
  ggplot(aes(bin, signal, color = decile)) +
  geom_line(size = 0.6) +
  scale_color_manual(values = brewer.pal(9, "Blues")[3:9]) +
  ylab(expression(Log[2]~"ChIP / Control")) +
  scale_y_continuous(limits = c(-0.4, 1.2), breaks = c(-0.4, 0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(breaks = c(0, 200, 400, 600), labels = c(2, "TSS", "TTS", 2)) +
  facet_wrap(~chip) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11)
  )
fig

ggsave(fig, filename = "figures/metaplots/expression_legend.pdf",
       height = 5, width = 14, units = "cm", dpi = 1200)
