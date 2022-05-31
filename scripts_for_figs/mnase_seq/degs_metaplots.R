library(dplyr)
library(stringr)
library(tidyr)


setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

data <- read.table("data/direct_targets_shift_tss.tab", skip = 2)
data[1:6, 1:6]

colnames(data) <- c("background", "genes", c(1:400))


unique(data$background)
data2 <- data %>%
  pivot_longer(cols = as.character(c(1:400)), names_to = "bin", values_to = "value") %>%
  mutate(bin = as.integer(bin)) %>%
  filter(!(background %in% c("tpf-col0.pois_diff", "chr-col0.pois_diff"))) %>%
  mutate(reg = if_else(genes %in% c("direct_up_shift_tss.bed", "direct_up_noshift_tss.bed"),
                       "up", if_else(genes %in% c("direct_down_shift_tss.bed", "direct_down_noshift_tss.bed"), 
                                     "down", "n.s."))) %>%
  mutate(shift = if_else(genes %in% c("direct_up_shift_tss.bed", "direct_down_shift_tss.bed", "direct_shift_nodeg_tss.bed"), T, F))

table(data2$background)

data2 %>%
  mutate(shift = factor(shift, levels = c(T, F))) %>%
  mutate(reg = factor(reg, levels = c("up", "down", "n.s."))) %>%
  mutate(wt = if_else(background == "col0.Fnor.smooth", T, F)) %>%
  #filter(shift == T) %>%
  mutate(background = factor(background, levels = c("col0.Fnor.smooth", "tpf.Fnor.smooth", "chr.Fnor.smooth"))) %>%
  ggplot(aes(bin, value)) +
  geom_line(aes(color = background)) +
  facet_wrap(reg~shift, ncol = 2, scales = "free_x") +
  ylab("Norm. read count") +
  xlab("Distance to TSS (kb)") +
  scale_x_continuous(breaks = c(0, 200, 400), labels = c("-0.2", "TSS", "0.2")) +
  scale_y_continuous(breaks = c(1000, 3000, 5000)) +
  scale_color_manual(values = c("grey", brewer.pal(3, "Dark2")[3:2])) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))

ggsave("figures/degs_mnase/test.pdf", height = 9, width = 8, units = "cm", dpi = 1200)
