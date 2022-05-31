#!/usr/bin/env Rscript

# TPF target metaplots
# Jaime Perez Alemany
# 07/05/2022

setwd("~/../Dropbox/Bioinformatics/ChIP-seq TPF1 TPF2 CHR23/")

library(dplyr)
library(tidyr)


profile <- read.table("data/targets_matrix_2000.inter.tab", skip = 2)[, c(3:602)] %>%
  t() %>% data.frame()

colnames(profile) <- c("TPF1", "TPF2", "CHR23")  
profile$bin <- c(1:nrow(profile))

data <- profile %>%
  pivot_longer(cols = c(TPF1, TPF2, CHR23), names_to = "chip", values_to = "value") %>%
  mutate(chip = factor(chip, levels = c("TPF1", "TPF2", "CHR23")))

plot <- ggplot(data, aes(bin, value)) +
  geom_line(size = 1, aes(color = chip)) +
  ylab(expression(Log[2]~"ChIP / Control")) +
  scale_y_continuous(limits = c(-0.02667534, 1.2)) +
  scale_x_continuous(breaks = c(0, 200, 400, 600), labels = c("-2 kb", "TSS", "TTS", "2 kb")) +
  theme_bw() +
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(1, 3, 2)]) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11)
        )
plot
dir.create("figures/metaplots")
ggsave(plot, filename = "figures/metaplots/metaplots_template.pdf", 
       height = 6, width = 9, units = "cm", dpi = 1200)
