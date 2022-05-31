library(RColorBrewer)
library(eulerr)
library(ggtext)
library(ggplot2)
library(dplyr)

setwd("~/../Dropbox/Bioinformatics/RNA-seq tpf chr/")
tpf_df <- read.csv("deseq2/tpf_vs_col0.counts.csv")
chr_df <- read.csv("deseq2/chr_vs_col0.counts.csv")

# Volcano plots

blue <- brewer.pal(4, "Paired")[2]
red <- brewer.pal(4, "Paired")[6]
display.brewer.pal(6, "Paired")
tpf_df$mutant <- "tpf1 tpf2"
chr_df$mutant <- "chr12 chr23"

rbind

fc_lim <- max(abs(tpf_df$log2FoldChange))
fc_lim
max(-log10(tpf_df$padj %>% na.omit()))

table(chr_df$deg)

template_tpf <- tpf_df %>%
  na.omit() %>%
  ggplot(aes(log2FoldChange, -log10(padj))) +
  theme_bw() +
  scale_color_manual(values = c(4, "grey", 2)) +
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, length.out = 5)) +
  scale_x_continuous(limits = c(-10.8, 10.8)) +
  ylab(expression(-Log[10]~italic(P))) +
  xlab("Log<sub>2</sub> (*tpf* / Col-0)") +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_markdown(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))


template_tpf
plot_tpf <- template_tpf +
  geom_point(aes(color = deg), alpha = 0.25, size = 2)

plot_tpf
template_tpf <- template_tpf +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
             color = "grey30", linetype = "dashed", size = 0.25) + 
  geom_hline(yintercept = -log10(0.05), 
             color = "grey30", linetype = "dashed", size = 0.25)

template_tpf
ggsave(template_tpf, filename = "figures/volcano_plot/tpf_template.pdf", 
       height = 5, width = 6, units = "cm", dpi = 1200)
ggsave(plot_tpf, filename = "figures/volcano_plot/tpf_plot.tiff",
       height = 5, width = 6, units = "cm", dpi = 1200)


template_chr <- chr_df %>% 
  na.omit() %>%
  ggplot(aes(log2FoldChange, -log10(padj))) +
  theme_bw() +
  scale_color_manual(values = c(4, "grey", 2)) +
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, length.out = 5)) +
  scale_x_continuous(limits = c(-10.8, 10.8)) +
  ylab(expression(-Log[10]~italic(P))) +
  xlab("Log<sub>2</sub> (*minu* / Col-0)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.title.x = element_markdown(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11))

plot_chr <- template_chr +
  geom_point(aes(color = deg), alpha = 0.25, size = 0.5)

template_chr <- template_chr +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
             color = "grey30", linetype = "dashed", size = 0.25) + 
  geom_hline(yintercept = -log10(0.05), 
             color = "grey30", linetype = "dashed", size = 0.25)

plot_chr
ggsave(template_chr, filename = "figures/volcano_plot/chr_template.pdf", 
       height = 5, width = 6, units = "cm", dpi = 1200)
ggsave(plot_chr, filename = "figures/volcano_plot/chr_plot.tiff",
       height = 5, width = 6, units = "cm", dpi = 1200)

dir.create("figures/venn_diagrams/")
pdf("figures/venn_diagrams/overlap_down.pdf")
venn(
  combinations = list(
    tpf = tpf_df[tpf_df$deg == "Down", "TAIR"] %>% na.omit(),
    chr = chr_df[chr_df$deg == "Down", "TAIR"] %>% na.omit()
  )
) %>% plot(edges = list(col = c("#1f78b4", "#a6cee3"), lex = 1.5), fills = NA, labels = NA)
dev.off()

pdf("figures/venn_diagrams/overlap_up.pdf")
venn(
  combinations = list(
    tpf = tpf_df[tpf_df$deg == "Up", "TAIR"] %>% na.omit(),
    chr = chr_df[chr_df$deg == "Up", "TAIR"] %>% na.omit()
  )
) %>% plot(edges = list(col = c("#e31a1c", "#fb9a99"), lex = 1.5), fills = NA, labels = NA)
dev.off()

colnames(tpf_df)[2:7] <- paste("tpf", colnames(tpf_df)[2:7], sep = "_")
colnames(chr_df)[2:7] <- paste("chr", colnames(chr_df)[2:7], sep = "_")

head(tpf_df)
chr_tpf <- merge(tpf_df, chr_df, by = "TAIR", all = T)

write.csv(chr_tpf, file = "deseq2/chr_tpf.integrative_dea.counts.csv")

pdf("figures/supp/correlation_fc.pdf", height = 4, width = 4)
smoothScatter(chr_tpf$tpf_log2FoldChange, chr_tpf$chr_log2FoldChange, xlim = c(-11, 11), ylim = c(-11, 11))
dev.off()

cor(chr_tpf$chr_log2FoldChange, chr_tpf$tpf_log2FoldChange)
chr_tpf %>%
  #filter(chr_deg != "NS" | tpf_deg != "NS") %>%
  ggplot(aes(tpf_log2FoldChange, chr_log2FoldChange)) +
  geom_point(size = 0.5, alpha = 0.25, color = brewer.pal(9, "Blues")[9]) +
  scale_x_continuous(limits = c(-10.8, 10.8)) +
  scale_y_continuous(limits = c(-10.8, 10.8)) +
  ylab("Log<sub>2</sub> (*minu1 minu2* / WT)") +
  xlab("Log<sub>2</sub> (*tpf1 tpf2* / WT)") +
  theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12))

