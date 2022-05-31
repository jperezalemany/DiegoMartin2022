#!/usr/bin/env Rscript

# Correlation between ChIP-seq replicates
setwd("~/../Dropbox/Bioinformatics/ChIP-seq TPF1 TPF2 CHR23/")

library(RColorBrewer)
library(scales)
library(cowplot)
library(dplyr)
library(ggplot2)

template <- ggplot() +
  theme_bw() +
  scale_x_continuous(
    limits = c(1.5, 10.5), expand = c(0.01, 0.01), breaks = seq(2, 10, by = 2)
  ) +
  scale_y_continuous(
    limits = c(1.5, 10.5), expand = c(0.01, 0.01), breaks = seq(2, 10, by = 2)
  ) +
  xlab("Line 1") +
  ylab("Line 2") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.4),
    axis.ticks = element_line(color = "black", size = 0.4),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11)
  )

template
tpf1 <- read.csv("counts/TPF1.csv") %>%
  select(contains("cpm"))

cor(tpf1$TPF1.4_cpm, tpf1$TPF1.11_cpm)


tpf1_plot <- template +
  geom_point(
    data = log2(tpf1), aes(TPF1.4_cpm, TPF1.11_cpm), size = 0.25, alpha = 0.25,
    color = brewer.pal(9, "Blues")[9]
  )

tpf1_plot
tpf2 <- read.csv("counts/TPF2.csv") %>%
  select(contains("cpm"))

cor(tpf2$TPF2.15_cpm, tpf2$TPF2.9_cpm)

tpf2_plot <- template +
  geom_point(
    data = log2(tpf2), aes(TPF2.9_cpm, TPF2.15_cpm), size = 0.25, alpha = 0.25,
    color = brewer.pal(9, "Blues")[9]
  )

pdf("figures/correlation/why_not.pdf", height = 4, width = 4)
smoothScatter(log2(tpf1), xlim = c(1.5, 10.5), ylim = c(1.5, 10.5), xlab = "Line 1", ylab = "Line 2",nbin = 100 )
dev.off()
pdf("figures/correlation/chr23.pdf", height = 4, width = 4)
smoothScatter(log2(chr23), xlim = c(1.5, 10.5), ylim = c(1.5, 10.5), xlab = "Line 1", ylab = "Line 2", nbin = 100)
dev.off()
pdf("figures/correlation/tpf2.pdf", height = 4, width = 4)
smoothScatter(log2(tpf2), xlim = c(1.5, 10.5), ylim = c(1.5, 10.5), xlab = "Line 1", ylab = "Line 2", nbin = 100)
dev.off()

chr23 <- read.csv("counts/CHR23_2.csv") %>%
  select(contains("cpm"))

cor(chr23$CHR23.7_2_cpm, chr23$CHR23.10_2_cpm)
head(chr23)

chr23_plot
chr23_plot <- template +
  geom_point(
    data = log2(chr23), aes(CHR23.7_2_cpm, CHR23.10_2_cpm),
    size = 0.25, alpha = 0.25, color = brewer.pal(9, "Blues")[9]
  )
chr23_plot
ggsave(chr23_plot, filename = "test.pdf", height = 5, width = 10, units = "cm", dpi = 1200)

plot <- plot_grid(tpf1_plot, tpf2_plot, chr23_plot, ncol = 3)
plot
dir.create("figures/correlation/", showWarnings = F)
ggsave(
  plot, filename = "figures/correlation/correlation_plots.pdf",
  height = 6, width = 18, units = "cm", dpi = 1200
)
ggsave(
  plot, filename = "figures/correlation/correlation_plots.tiff",
  height = 6, width = 18, units = "cm", dpi = 1200
)

plot_grid(tpf1_plot, tpf2_plot, chr23_plot, ncol = 3)
ggsave(filename = "figures/correlation/correlation_plots2.pdf", height = 6, width = 18, units = "cm", dpi = 1200)

template_labs <- template
empty <- plot_grid(
  template_labs, template_labs, template_labs, ncol = 3
)
ggsave(
  empty, filename = "figures/correlation/correlation_plots_empty.pdf",
  height = 6, width = 18, units = "cm", dpi = 1200
)
