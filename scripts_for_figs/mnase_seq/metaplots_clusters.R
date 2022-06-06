library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)


setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

data <- read.table("data/mnase_rna_chip_targets_notargets.tss_signal.tab", skip = 2)

colnames(data) <- c("sample", "cluster", c(1:400))

data[1:6, 1:6]
data2 <- data %>%
  pivot_longer(cols = as.character(c(1:400)), names_to = "bin", values_to = "signal") %>%
  mutate(bin = as.integer(bin)) %>%
  mutate(omic = str_remove(sample, "_(col0|tpf|chr)")) %>%
  mutate(background = str_remove(sample, "(mnase|rna)_")) %>%
  mutate(background = factor(background, levels = rev(c("col0", "tpf", "chr")))) %>%
  mutate(cluster = factor(cluster, levels = c("targets", "notargets")))

max(data2[data2$omic == "mnase", "signal"])
mnase <- data2 %>%
  filter(omic == "mnase") %>%
  mutate(background = factor(background, levels = c("col0", "tpf", "chr"))) %>%
  ggplot(aes(bin, signal, color = background)) +
  scale_x_continuous(breaks = c(0, 200, 400), labels = c(-0.2, "TSS", 0.2)) +
  scale_y_continuous(limits = c(1000, 5000), breaks = c(1000, 3000, 5000)) +
  xlab("Distance from TSS (Kb)") +
  geom_line(size = 0.75) +
  theme_bw() +
  facet_wrap(~cluster) +
  scale_color_manual(values = c("grey", brewer.pal(3, "Dark2")[3:2])) +
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

rna <- data2 %>%
filter(omic == "rna") %>%
  mutate(background = factor(background, levels = c("col0", "tpf", "chr"))) %>%
  ggplot(aes(bin, signal, color = background)) +
  scale_x_continuous(breaks = c(0, 200, 400), labels = c(-0.2, "TSS", 0.2)) +
  scale_y_continuous(limits = c(0, 5), breaks = c(0, 2.5, 5)) +
  xlab("Distance from TSS (Kb)") +
  geom_line(size = 0.75) +
  theme_bw() +
  facet_wrap(~cluster) +
  scale_color_manual(values = c("grey", brewer.pal(3, "Dark2")[3:2])) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))


plot_grid(mnase, rna, ncol = 1, align = "v")
dir.create("figures/clusters")
ggsave("figures/clusters/metaplot_clusters.pdf", height = 6, width = 9, units = "cm", dpi = 1200)
