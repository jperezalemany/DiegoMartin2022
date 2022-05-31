
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

data <- read.table("data/plus_one_byreps.tab", skip = 2)

colnames(data) <- c("group", "genes", c(1:150))
data$background <- factor(
  c("col0", "col0", "tpf", "tpf", "chr", "chr", "choi"), levels = c("col0", "tpf", "chr", "choi")
)
data$rep <- c(as.character(c(1, 2, 1, 2, 1, 2, 1)))

data

data %>%
  pivot_longer(as.character(c(1:150)), values_to = "value", names_to = "bin") %>%
  mutate(bin = as.integer(bin)) %>%
  ggplot(aes(bin, value, color = rep)) +
  facet_wrap(~background, nrow = 1) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[2:3]) +
  scale_x_continuous(labels = c(-0.5, 0, 0.5, 1)) +
  scale_y_continuous(limits = c(1500, 6500), breaks = c(1500, 4000, 6500)) +
  xlab("Distance to +1 dyad (kb)") +
  ylab("Norm. counts") +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))

ggsave(filename = "figures/supp/phasing.pdf", height = 5, width = 20, units = "cm", dpi = 1200)

