library(clusterProfiler)
library(org.At.tair.db)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)

setwd("~/../Dropbox/Bioinformatics/RNA-seq tpf chr")

tpf_df <- read.csv("deseq2/tpf_vs_col0.counts.csv")
chr_df <- read.csv("deseq2/chr_vs_col0.counts.csv")

table(tpf_df$deg)

tpf_ora <- enrichGO(
  gene = tpf_df[tpf_df$deg != "NS", "TAIR"],
  OrgDb = org.At.tair.db,
  ont = "BP", keyType = "TAIR",
  pvalueCutoff = 1, qvalueCutoff = 1
)

chr_ora <- enrichGO(
  gene = chr_df[chr_df$deg != "NS", "TAIR"],
  OrgDb = org.At.tair.db,
  ont = "BP", keyType = "TAIR",
  pvalueCutoff = 1, qvalueCutoff = 1
)

tpf_ora_sig <- tpf_ora %>% filter(p.adjust < 0.05) %>% simplify(cutoff = 0.5)
chr_ora_sig <- chr_ora %>% filter(p.adjust < 0.05) %>% simplify(cutoff = 0.5)


chr_cut <- chr_ora_sig@result %>%
  dplyr::distinct(p.adjust, .keep_all = T) %>%
  arrange(p.adjust) %>%
  head(n = 5) %>%
  tail(n = 1)

tpf_cut <- chr_ora_sig@result %>%
  dplyr::distinct(p.adjust, .keep_all = T) %>%
  arrange(p.adjust) %>%
  head(n = 5) %>%
  tail(n = 1)


chr_cut
pathwats <- c(
  chr_ora_sig@result[tpf_ora_sig@result$p.adjust <= chr_cut$p.adjust, "Description"],
  tpf_ora_sig@result[tpf_ora_sig@result$p.adjust <= tpf_cut$p.adjust, "Description"]
) %>% unique()

pathwats
data <- rbind(
  chr_ora@result %>%
  filter(Description %in% pathwats) %>%
  mutate(mutant = "chr12_chr23"),
  tpf_ora@result %>%
  filter(Description %in% pathwats) %>%
  mutate(mutant = "tpf1_tpf2")
) %>%
  dplyr::select(mutant, Description, p.adjust) %>%
  dplyr::rename(padj = p.adjust)

levels <- data %>%
  mutate(padj = -log10(padj)) %>%
  pivot_wider(id_cols = Description, names_from = mutant, values_from = padj)

levels <- levels %>%
  arrange(desc((tpf1_tpf2 + chr12_chr23)/2))

head(levels)
colnames(chr_ora@result)
library(ComplexHeatmap)

data %>%
  mutate(Description = factor(Description, levels = rev(levels$Description))) %>%
  #mutate(mutant = factor(mutant, levels = c("tpf1_tpf2", "chr12_chr23"))) %>%
  ggplot(aes(Description, -log10(padj), fill = mutant)) +
  scale_fill_manual(values = brewer.pal(3, "Set2")[2:3]) +
  geom_bar(stat = "identity", position = "dodge", width = 0.75) +
  #geom_hline(yintercept = -log10(0.05)) +
  scale_y_continuous(expand = c(0, 0, 0.05, 0.05), limits = c(0, 25)) +
  coord_flip() +
  ylab(expression(-Log[10]~italic(P))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11)
  )

ggsave(filename = "figures/supp/go_barplot.pdf", height = 15, width = 15, units = "cm", dpi = 1200)

library(circlize)

col_fun <- colorRamp2(breaks = c(0, 30), colors = c("white", "blue4"))

pdf("figures/supp/go_heatmap.pdf", height = 10, width = 8)
data %>%
  #dplyr::rename(padj = p.adjust) %>%
  mutate(padj = -log10(padj)) %>%
  pivot_wider(id_cols = Description, names_from = mutant, values_from = padj) %>%
  tibble::column_to_rownames(var = "Description") %>%
  as.matrix() %>%
  Heatmap(col = col_fun)
dev.off()

data %>%
  ggplot(aes(as.numeric(-log10(p.adjust)), Description)) +
  geom_bar(stat = "identity")


data %>%
  ggplot(aes())
