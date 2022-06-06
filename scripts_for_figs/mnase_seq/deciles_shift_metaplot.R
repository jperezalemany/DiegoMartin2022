setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)


data <- read.table("data/deciles_shift_plus.tab", skip = 2)
data[1:6, 1:6]
colnames(data) <- c("sample", "cluster", c(1:1000))

data2 <- data %>%
  pivot_longer(cols = as.character(c(1:1000)), names_to = "bin", values_to = "signal") %>%
  mutate(bin = as.integer(bin)) %>%
  mutate(omic = str_remove(sample, "_(col0|tpf|chr)")) %>%
  mutate(background = str_remove(sample, "(mnase|rna)_")) %>%
  mutate(background = factor(background, levels = c("tpf", "chr"))) %>%
  mutate(cluster = factor(cluster, levels = c(1:5)))

plot1 <- data2 %>%
  filter(omic == "rna") %>%
  ggplot(aes(bin, signal, color = cluster)) +
  facet_wrap(~background, scales = "free_y") +
  geom_line(size = 0.75) +
  scale_color_manual(values = rev(brewer.pal(9, "Blues")[4:8])) +
  scale_x_continuous(breaks = c(0, 500, 1000), labels = c(-0.5, 0, 0.5)) +
  scale_y_continuous(breaks = c(0, 0.6, 1.2, 1.8), limits = c(-0.1, 1.8)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))

min(data2[data2$omic == "mnase", "signal"])
plot2 <- data2 %>%
  filter(omic == "mnase") %>%
  ggplot(aes(bin, signal, color = cluster)) +
  facet_wrap(~background, scales = "free_y") +
  geom_line(size = 0.75) +
  scale_color_manual(values = rev(brewer.pal(9, "Blues")[4:8])) +
  scale_y_continuous(breaks = c(-250, 0, 250, 500), limits = c(-290, 500)) +
  scale_x_continuous(breaks = c(0, 500, 1000), labels = c(-0.5, 0, 0.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "right",
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))
plot2

ggsave(plot2,
       filename = "figures/metaplots_shift_targets/legend.pdf",
       height = 6, width = 12, units = "cm", dpi = 1200)
library(cowplot)

plot_grid(plot2, plot1, ncol = 1, align = "v")
ggsave("figures/metaplots_shift_targets/deciles.pdf", height = 8, width = 12, units = "cm", dpi = 1200)
