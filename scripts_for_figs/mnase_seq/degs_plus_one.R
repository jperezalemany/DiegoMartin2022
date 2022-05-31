
library(dplyr)

setwd("~/../Dropbox/Bioinformatics/")

targets <- read.table("ChIP-seq TPF1 TPF2 CHR23/bedfiles/TPF.inter_targets.bed")
indirect_targets <- read.table("ChIP-seq TPF1 TPF2 CHR23/bedfiles/TPF.inter_notargets.bed")

degs <- read.csv("RNA-seq tpf chr/deseq2/chr_tpf.integrative_dea.counts.csv", row.names = 1) %>%
  mutate(target = if_else(TAIR %in% targets$V4, "target", "NS")) %>%
  mutate(target = if_else(TAIR %in% indirect_targets$V4, "notarget", target))



plus_one <- read.table("MNase-seq tpf chr/bedfiles/col0_plus_one.bed")

direct_down <- degs %>%
  filter(chr_deg == "Down" & tpf_deg == "Down" & target == "target")

direct_up <- degs %>%
  filter(tpf_deg == "Up" & chr_deg == "Up" & target == "target")

indirect_down <- degs %>%
  filter(chr_deg == "Down" & tpf_deg == "Down" & target == "notarget")

indirect_up <- degs %>%
  filter(chr_deg == "Up" & tpf_deg == "Up" & target == "notarget")

direct_static <- targets %>%
  filter(!(V4 %in% c(direct_down$TAIR, direct_up$TAIR)))

  
gene_bed <- read.table("Araport11/TSS_adjusted.bed")

dir.create("MNase-seq tpf chr/plus_one/")
plus_one %>%
  filter(V11 %in% direct_up$TAIR) %>%
  write.table("MNase-seq tpf chr/plus_one/targets_vs_degs/direct_up.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)

gene_bed %>%
  filter(V4 %in% direct_up$TAIR) %>%
  write.table("MNase-seq tpf chr/plus_one/targets_vs_degs/direct_up_tss.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)

plus_one %>%
  filter(V11 %in% direct_down$TAIR) %>%
  write.table("MNase-seq tpf chr/plus_one/targets_vs_degs/direct_down.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)

gene_bed %>%
  filter(V4 %in% direct_down$TAIR) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_down_tss.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)

plus_one %>%
  filter(V11 %in% direct_static$TAIR) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_static.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)

gene_bed %>%
  filter(V4 %in% direct_static$V4) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_static_tss.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)

#shift <- read.csv("MNase-seq tpf chr/processed_tables/chr_tpf.integrated_shift.csv")
plus_one_ann <- read.csv("MNase-seq tpf chr/processed_tables/plus_one_data.csv") %>%
  mutate(mean_shift = (chr2control_dis + tpf2control_dis)/2) %>%
  filter(col0_fuzz_pval < log10(0.05), tpf_fuzz_pval < log10(0.05), chr_fuzz_pval < log10(0.05)) %>%
  filter(tpf_smt_pval < log10(0.05), chr_smt_pval < log10(0.05))

shifted_plus_one <- plus_one_ann %>%
  filter(mean_shift < -5) %>%
  filter(chr2control_dis < -5 & tpf2control_dis < -5) %>%
  filter(target == T)

eulerr::euler(
  list(
    up = direct_up$TAIR,
    down = direct_down$TAIR,
    shift = shifted_plus_one$gene_name,
    targets = targets$V4
  )
) %>% plot(quantities = T)


# Final categories
down_shift <- intersect(direct_down$TAIR, shifted_plus_one$gene_name)
up_shift <- intersect(direct_up$TAIR, shifted_plus_one$gene_name)

down_noshift <- setdiff(direct_down$TAIR, shifted_plus_one$gene_name)
up_noshift <- setdiff(direct_up$TAIR, shifted_plus_one$gene_name)

shift <- setdiff(shifted_plus_one$gene_name, c(down_shift, up_shift))
noshift <- setdiff(targets$V4, c(shift, up_shift, down_shift, down_noshift, up_noshift))

data.frame(
  deg = c("Direct up", "Direct up","Direct down", "Direct down", "Direct NS", "Direct NS"),
  change = c("upstream", "other", "upstream", "other", "upstream", "other"),
  number = c(
    length(up_shift), length(up_noshift), length(down_shift), length(down_noshift),
    length(shift), length(noshift)
  ),
  from = c(
    nrow(direct_up), nrow(direct_up), nrow(direct_down), nrow(direct_down),
    length(shift) + length(noshift), length(shift) + length(noshift)
  )
) %>%
  mutate(percent = number / from * 100) %>%
  mutate(deg = factor(deg, levels = c("Direct NS", "Direct down", "Direct up"))) %>%
  ggplot(aes(deg, percent, fill = change)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ylab("Upstream shift at +1 dyad (%)") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        axis.line = element_line(color = "black", size = 0.4),
        axis.ticks = element_line(color = "black", size = 0.4),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 11)) +
  scale_fill_manual(values = brewer.pal(9, "Paired"))

dir.create("MNase-seq tpf chr/figures/degs_mnase/")
ggsave(filename = "MNase-seq tpf chr/figures/degs_mnase/percents.pdf",
       height = 4, width = 10, dpi = 1200, units = "cm")

gene_bed %>%
  filter(V4 %in% down_shift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_down_shift_tss.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)
gene_bed %>%
  filter(V4 %in% down_noshift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_down_noshift_tss.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)

plus_one %>%
  filter(V11 %in% down_shift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_down_shift.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)
plus_one %>%
  filter(V11 %in% down_noshift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_down_noshift.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)

gene_bed %>%
  filter(V4 %in% up_shift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_up_shift_tss.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)
gene_bed %>%
  filter(V4 %in% up_noshift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_up_noshift_tss.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)

plus_one %>%
  filter(V11 %in% up_shift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_up_shift.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)
plus_one %>%
  filter(V11 %in% up_noshift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_up_noshift.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)


plus_one %>%
  filter(V11 %in% shift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_shift_nodeg.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)
gene_bed %>%
  filter(V4 %in% shift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_shift_nodeg_tss.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)

plus_one %>%
  filter(V11 %in% noshift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_noshift_nodeg.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)
gene_bed %>%
  filter(V4 %in% noshift) %>%
  write.table("MNase-seq tpf chr/plus_one/direct_targets_vs_shift/direct_noshift_nodeg_tss.bed",
              sep = "\t", quote = F, row.names = F, col.names = F)


