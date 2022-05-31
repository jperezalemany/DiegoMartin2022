
library(dplyr)
library(eulerr)
library(ggplot2)
library(tidyr)
library(DESeq2)
library(RColorBrewer)

setwd("~/../Dropbox/Bioinformatics/RNA-seq tpf chr//")

tpf <- read.csv("deseq2/tpf_vs_col0.counts.csv")

tpf %>%
  write.table("deseq2/tpf_vs_col0.counts.tab", quote = F, row.names = F, sep = "\t")

read.csv("deseq2/chr_vs_col0.counts.csv") %>%
  write.table("deseq2/chr_vs_col0.counts.tab", sep = "\t", row.names = F, quote = F)

degs <- read.csv("../RNA-seq tpf chr/deseq2/chr_tpf.integrative_dea.counts.csv")
targets <- read.table("../ChIP-seq TPF1 TPF2 CHR23/bedfiles/TPF.inter_targets.bed")
no_targets <- read.table("../ChIP-seq TPF1 TPF2 CHR23/bedfiles/TPF.inter_notargets.bed")

table <- degs %>%
  mutate(target = if_else(TAIR %in% targets$V4, "target",
                          if_else(TAIR %in% no_targets$V4, "not_target", "discarded"))) %>%
  filter(!is.na(chr_padj)) %>%
  filter(!is.na(tpf_padj))

table[table$chr_deg == "Down" & table$tpf_deg == "Down", "deg"] <- "common_down"
table[table$chr_deg == "Up" & table$tpf_deg == "Up", "deg"] <- "common_up"
table[table$chr_deg == "Down" & table$tpf_deg == "NS", "deg"] <- "chr_down"
table[table$chr_deg == "Up" & table$tpf_deg == "NS", "deg"] <- "chr_up"
table[table$tpf_deg == "Down" & table$chr_deg == "NS", "deg"] <- "tpf_down"
table[table$tpf_deg == "Up" & table$chr_deg == "NS", "deg"] <- "tpf_up"
table[table$chr_deg == "Up" & table$tpf_deg == "Down", "deg"] <- "diff"
table[table$chr_deg == "Down" & table$tpf_deg == "Up", "deg"] <- "diff"
table[table$chr_deg == "NS" & table$tpf_deg == "NS", "deg"] <- "common_not_deg"

table[table$chr_deg == "Down" & table$tpf_deg == "Down", "deg2"] <- "common"
table[table$chr_deg == "Up" & table$tpf_deg == "Up", "deg2"] <- "common"
table[table$chr_deg == "Down" & table$tpf_deg == "NS", "deg2"] <- "chr"
table[table$chr_deg == "Up" & table$tpf_deg == "NS", "deg2"] <- "chr"
table[table$tpf_deg == "Down" & table$chr_deg == "NS", "deg2"] <- "tpf"
table[table$tpf_deg == "Up" & table$chr_deg == "NS", "deg2"] <- "tpf"
table[table$chr_deg == "Up" & table$tpf_deg == "Down", "deg2"] <- "diff"
table[table$chr_deg == "Down" & table$tpf_deg == "Up", "deg2"] <- "diff"

table(table$deg, table$target)

dds <- readRDS("../RNA-seq tpf chr/deseq2/dds.RDS")
mcols(dds)$basepairs <- mcols(dds)$basepair

exp <- fpkm(dds) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "TAIR") %>%
  mutate(mean_col0 = (col0_rep1 + col0_rep2 + col0_rep3) / 3) %>%
  mutate(mean_chr = (chr_rep1 + chr_rep2 + chr_rep3) / 3) %>%
  mutate(mean_tpf = (tpf_rep1 + tpf_rep2 + tpf_rep3) / 3) %>%
  select(TAIR, contains("mean")) %>%
  pivot_longer(contains("mean"), values_to = "fpkm", names_to = "mutant")


merge(table, exp, by = "TAIR") %>%
  pivot_longer(contains("log2FoldChange"), names_to = "mutant2", values_to = "FC") %>%
  filter(deg == "common_down" | deg == "common_up") %>%
  mutate(deg = if_else(deg == "common_down", "Downregulated", "Upregulated")) %>%
  mutate(deg = factor(deg, levels = c("Upregulated", "Downregulated"))) %>%
  mutate(mutant = if_else(mutant == "mean_col0", "wild type", if_else(mutant == "mean_tpf", "tpf1 tpf2", "minu1 minu2"))) %>%
  filter(target == "target") %>%
  mutate(mutant = factor(mutant, levels = c("wild type", "tpf1 tpf2", "minu1 minu2"))) %>%
  ggplot(aes(deg, log2(fpkm))) +
  stat_boxplot(geom = "errorbar", aes(color = mutant), size = 0.5) +
  scale_color_manual(values = brewer.pal(3, "Set2")[c(1, 3, 2)]) +
  geom_boxplot(aes(color = mutant), outlier.shape = 1, outlier.size = 1, size = 0.5) +
  ylab(expression(Log[2]~FPKM)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12))

ggsave(filename = "figures/test_boxplot.pdf", height = 6, width = 9, units = "cm", dpi = 1200)

genes <- read.table("../Araport11/all_genes.bed")

dir.create("bedfiles")
genes %>%
  filter(V4 %in% table[table$deg == "common_down" & table$target == "target", "TAIR"]) %>%
  write.table("bedfiles/common_down.bed", sep = "\t", row.names = F, col.names = F,
              quote = F)

genes %>%
  filter(V4 %in% table[table$deg == "common_up" & table$target == "target", "TAIR"]) %>%
  write.table("bedfiles/common_up.bed", sep = "\t", row.names = F, col.names = F,
              quote = F)

genes %>%
  filter(V4 %in% table[table$deg2 == "common" & table$target == "target", "TAIR"]) %>%
  write.table("bedfiles/common_all.bed", sep = "\t", row.names = F, col.names = F,
              quote = F)

euler(
  list(
    down = table[table$deg == "common_down", "TAIR"] %>% unique(),
    up = table[table$deg == "common_up", "TAIR"] %>% unique(),
    targets = targets$V4
  )
) %>% plot(quantities = T)

data.frame(
  deg = c("Down", "Down", "Up", "Up"),
  target = c("TPF-occ.", "TPF-free", "TPF-occ.", "TPF-free"),
  n = c(531, 1305, 427, 354)
) %>%
  ggplot(aes(target, n, fill = deg)) +
  #scale_fill_manual(values = c("red3", "blue3")) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_bw()
  
  



