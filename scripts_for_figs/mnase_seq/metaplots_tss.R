library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")


data <- read.table("data/mnase_rna_chip_targets_notargets.tss.tab", skip = 2)
colnames(data) <- c("sample", "cluster", c(1:150))
data[1:6, 1:6]

data2 <- data %>%
  pivot_longer(cols = as.character(c(1:150)), names_to = "bin", values_to = "signal") %>%
  mutate(bin = as.integer(bin)) %>%
  mutate(omic = str_remove(sample, "_(col0|tpf|chr)")) %>%
  mutate(background = str_remove(sample, "(mnase|rna)_")) %>%
  mutate(background = factor(background, levels = rev(c("col0", "tpf", "chr")))) %>%
  mutate(cluster = factor(cluster, levels = c("targets", "notargets")))


tss_plot <- data2 %>%
  filter(omic == "rna") %>%
  mutate(background = factor(background, levels = c("tpf", "chr"))) %>%
  ggplot(aes(bin, signal, color = cluster)) +
  xlab("Distance from TSS (kb)") +
  geom_line(size = 0.75) +
  scale_y_continuous(breaks = c(-0.4, 0, 0.4, 0.8), limits = c(-0.4, 0.8)) +
  scale_x_continuous(breaks = c(0, 50, 100, 150), labels = c(-0.5, 0, 0.5, 1)) +
  theme_bw() +
  facet_wrap(~background, scales = "free", ncol = 1) +
  scale_color_manual(values = brewer.pal(3, "Dark2")) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))

ggsave(tss_plot, filename = "figures/metaplots_shift_targets/test_tss.pdf", height = 6.5, width = 4.5, units = "cm", dpi = 1200)
