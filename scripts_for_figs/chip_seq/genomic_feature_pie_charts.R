#!/usr/bin/env Rscript

# Genome annotation pie charts
# Jaime Perez Alemany
# 25/04/2022

setwd("~/../Dropbox/Bioinformatics/ChIP-seq TPF1 TPF2 CHR23/")

library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

get_feature_frequency <- function(annotated_peaks, column, name) {
  result <- annotated_peaks %>%
    group_by(!!as.name(column)) %>%
    summarise(n = n()) %>%
    rename(feature = !!as.name(column)) %>%
    mutate(freq = n / sum(n) * 100) %>%
    mutate(name = name)
  return(result)
}

tpf1 <- read.table("peaks/TPF1_summits.annotation.bed") %>%
  filter(V17 != "gene_inside_peak") %>%
  select(V4, V17) %>%
  distinct() %>%
  get_feature_frequency(column = "V17", name = "TPF1")

tpf2 <- read.table("peaks/TPF2_summits.annotation.bed") %>%
  filter(V17 != "gene_inside_peak") %>%
  select(V4, V17) %>%
  distinct() %>%
  get_feature_frequency(column = "V17", name = "TPF2")

chr23 <- read.table("peaks/CHR23_2_summits.annotation.bed") %>%
  filter(V17 != "gene_inside_peak") %>%
  select(V4, V17) %>%
  distinct() %>%
  get_feature_frequency(column = "V17", name = "MINU2")

genome <- read.table("bedfiles/Araport11.bed", sep = "\t") %>%
  select(V1, V2, V3, V4) %>%
  distinct() %>%
  mutate(width = V3 - V2) %>%
  group_by(V4) %>%
  summarise(n = sum(width)) %>%
  rename(feature = V4) %>%
  mutate(freq = n / sum(n) * 100) %>%
  mutate(name = "Genome")


features <- c(
  "Exon", "Intron", "5' UTR", "3' UTR", "Proximal promoter", "Terminator",
  "Distal promoter", "Intergenic"
)

groups <- c("TPF1", "TPF2", "MINU2", "Genome")

data <- rbind(tpf1, tpf2, chr23, genome) %>%
  mutate(name = factor(name, levels = groups)) %>%
  mutate(feature = as.character(feature)) %>%
  mutate(feature = str_to_title(feature)) %>%
  mutate(feature = str_replace_all(feature, "_", " ")) %>%
  mutate(feature = if_else(feature == "Five prime utr", "5' UTR", feature)) %>%
  mutate(feature = if_else(feature == "Three prime utr", "3' UTR", feature)) %>%
  mutate(feature = factor(feature, levels = features))

ggplot(data, aes(x = "", y = freq, fill = feature)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.4) +
  coord_polar("y", start= 0) +
  facet_wrap(~ name, ncol = 2) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  theme_void() +
  theme(legend.position = "right",
        plot.background = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 11),
        legend.title = element_blank())

dir.create("figures/genomic_features_pie_charts", showWarnings = F)
ggsave(
  filename = "figures/genomic_features_pie_charts/pie_charts_over_bases.pdf",
  height = 7, width = 12, dpi = 300, units = "cm"
)

data %>%
  write.csv("figures/genomic_features_pie_charts/pie_charts_over_bases.csv")
