#!/usr/bin/Rscript

# Metaplots MNase targets vs no targets

setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

library(dplyr)
library(stringr)
library(RColorBrewer)

data <- read.table("data/plus_one.wtmfuzz_targets_notargets.tab", sep = "\t", skip = 2)
colnames(data) <- c("background", "genes", paste0("bin_", c(1:600)))
data[1:6, 1:6]

data[data$background == "tpf1_tpf2", "background"] <- "*tpf*"
data[data$background == "minu1_minu2", "background"] <- "*minu*" 
data[data$genes == "targets", "genes"] <- "TPF-occ."
data[data$genes == "notargets", "genes"] <- "TPF-free" 

data <- data %>%
  pivot_longer(
    cols = contains("bin"), names_to = "bin", values_to = "signal"
  ) %>%
  mutate(bin = str_remove(bin, "bin_")) %>%
  mutate(bin = as.integer(bin)) %>%
  mutate(wt = if_else(background == "WT", T, F))

head(data)

max(data$signal)
plot_targets <- data %>%
  mutate(background = factor(background, levels = c("WT", "*tpf*", "*minu*"))) %>%
  filter(genes == "TPF-occ.") %>%
  ggplot(aes(bin, signal)) +
  geom_line(aes(color = background), size = 0.75) +
  scale_x_continuous(labels = c(-200, 0, 200, 400)) +
  scale_y_continuous(limits = c(1000, 7000), breaks = c(1000, 4000, 7000)) +
  xlab("Distance from +1 dyad (bp)") +
  ylab("Normalized read count") +
  scale_linetype_manual(values = c(1, 2)) +
  scale_color_manual(values = c("grey", brewer.pal(3, "Dark2")[3:2])) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        strip.text = element_markdown(),
        strip.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11)
  )
plot
library(ggtext)
ggsave(plot, filename = "figures/metaplots_shift_targets/test1.pdf", height = 5, width = 10, units = "cm", dpi = 1200)

plot_notargets <- data %>%
  filter(genes == "TPF-free") %>%
  ggplot(aes(bin, signal)) +
  geom_line(aes(color = background), size = 0.75) +
  scale_x_continuous(labels = c(-200, 0, 200, 400)) +
  scale_y_continuous(limits = c(1000, 7000), breaks = c(1000, 4000, 7000)) +
  xlab("Distance from +1 dyad (bp)") +
  ylab("Normalized read count") +
  scale_linetype_manual(values = c(1, 2)) +
  scale_color_manual(values = c(brewer.pal(3, "Dark2")[2:3], "grey")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11)
  )
plot
ggsave(plot, filename = "figures/metaplots_shift_targets/test2.pdf", height = 5, width = 10, units = "cm", dpi = 1200)


library(cowplot)
plot_notargets

plot_grid(plot_targets, plot_notargets, ncol = 1)
ggsave(filename = "figures/metaplots_shift_targets/targets_notargets.pdf", height = 8, width = 10, units = "cm", dpi = 1200)
