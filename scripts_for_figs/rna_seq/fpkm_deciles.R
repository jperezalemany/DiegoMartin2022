
library(DESeq2)
setwd("~/../Dropbox/Bioinformatics/RNA-seq tpf chr/")

dds_counts <- readRDS("deseq2/dds.RDS")
mcols(dds_counts)$basepairs <- mcols(dds_counts)$basepair

#dds_salmon <- readRDS("deseq2/dds_salmon.RDS")

library(dplyr)
exp_counts <- fpkm(dds_counts) %>% data.frame() %>%
  mutate(mean_col0 = (col0_rep1 + col0_rep2 + col0_rep3) / 3) %>%
  tibble::rownames_to_column(var = "TAIR") %>%
  select(TAIR, mean_col0)

#exp_salmon <- fpkm(dds_salmon) %>% data.frame() %>%
#  mutate(mean_col0 = (col0_rep1 + col0_rep2 + col0_rep3) / 3) %>%
#  tibble::rownames_to_column(var = "TAIR") %>%
#  select(TAIR, mean_col0)
 
greylist <- read.table("../ChIP-seq TPF1 TPF2 CHR23/bedfiles/greylist_ann.bed")

genes <- read.table("../Araport11/genes_mRNA.bed") %>%
  filter(V1 != "ChrC", V1 != "ChrM") %>%
  filter(!(V4 %in% greylist$V7))

# For salmon:
genes_salmon <- merge(genes, exp_salmon, by.x = "V4", by.y = "TAIR", all.x = T) %>%
  arrange(mean_col0)

not_exp_salmon <- genes_salmon[is.na(genes_salmon$mean_col0) | genes_salmon$mean_col0 == 0, ]
not_exp_salmon[, "decile"] <- 0

exp_salmon <- genes_salmon[!(genes_salmon$V4 %in% not_exp_salmon$V4), ]
decile_pos <- round(seq(1, nrow(exp_salmon), length.out = 6))

for (i in c(1:(length(decile_pos) -1))) {
  exp_salmon[decile_pos[i]:decile_pos[i+1], "decile"] <- i
}

genes_salmon <- rbind(not_exp_salmon, exp_salmon)
genes_salmon %>%
  mutate(log2_exp = log2(mean_col0)) %>% group_by(decile) %>%
  summarise(low = min(log2_exp), high = max(log2_exp))

table(genes_salmon$decile)

for (d in unique(genes_salmon$decile)) {
  genes_salmon %>%
    filter(decile == d) %>%
    select(V1, V2, V3, V4, V5, V6, mean_col0) %>%
    write.table(paste("fpkm_deciles/salmon_", d, ".bed", sep = ""),
                row.names = F, col.names = F, quote = F, sep = "\t")
}

# For gene-counts:
genes_counts <- merge(genes, exp_counts, by.x = "V4", by.y = "TAIR", all.x = T) %>%
  arrange(mean_col0)

not_exp_counts <- genes_counts[is.na(genes_counts$mean_col0) | genes_counts$mean_col0 == 0, ]
not_exp_counts[, "decile"] <- 0

exp_counts <- genes_counts[!(genes_counts$V4 %in% not_exp_counts$V4), ]
decile_pos <- round(seq(1, nrow(exp_counts), length.out = 6))

for (i in c(1:(length(decile_pos) -1))) {
  exp_counts[decile_pos[i]:decile_pos[i+1], "decile"] <- i
}

genes_counts <- rbind(not_exp_counts, exp_counts)
head(genes_counts)

genes_counts %>%
  mutate(log2exp = log2(mean_col0)) %>%
  mutate(decile = decile + 1) %>%
  arrange(desc(mean_col0)) %>%
  write.table("data/deciles_exp.bed", col.names = F, row.names = F, quote = F, sep = "\t")
  group_by(decile) %>%
  summarise(from = max(log2exp), to = min(log2exp))

for (d in unique(genes_counts$decile)) {
  genes_counts %>%
    filter(decile == d) %>%
    select(V1, V2, V3, V4, V5, V6, mean_col0) %>%
    write.table(paste("fpkm_deciles/counts_", d, ".bed", sep = ""),
                row.names = F, col.names = F, quote = F, sep = "\t")
}


