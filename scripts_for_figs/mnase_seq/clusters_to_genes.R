
library(dplyr)

bed_fields <- c("chr", "start", "end", "name", "score", "strand")

setwd("~/../Dropbox/Bioinformatics/MNase-seq tpf chr/")

clusters <- read.table("data/mnase_rna_chip_targets.tab") %>%
  select(V4, V13)
colnames(clusters) <- c("nucleosome", "cluster")

head(plus_one)

plus_one <- read.table("bedfiles/col0_plus_one.bed") %>%
  select(V4, V11)
colnames(plus_one) <- c("nucleosome", "gene")

genes_bed <- read.table("../Araport11/TSS_adjusted.bed") %>%
  select(V1, V2, V3, V4, V5, V6)
colnames(genes_bed) <- bed_fields

data <- merge(clusters, plus_one) %>% merge(genes_bed, by.x = "gene", by.y = "name")

dir.create("clusters")
data %>%
  filter(cluster == "cluster_1") %>%
  select(chr, start, end, gene, score, strand) %>%
  write.table("clusters/cluster_1_genes.bed", sep = "\t", quote = F, col.names = F, row.names = F)

data %>%
  filter(cluster == "cluster_2") %>%
  select(chr, start, end, gene, score, strand) %>%
  write.table("clusters/cluster_2_genes.bed", sep = "\t", quote = F, col.names = F, row.names = F)

data %>%
  filter(cluster == "cluster_3") %>%
  select(chr, start, end, gene, score, strand) %>%
  write.table("clusters/cluster_3_genes.bed", sep = "\t", quote = F, col.names = F, row.names = F)

data %>%
  filter(cluster == "cluster_3") %>%
  nrow()

chr23 %>%
  filter(V11 %in% data[data$cluster == "cluster_1", "gene"]) %>%
  select(V1, V2, V3, V4, V12, V13, V11) %>%
  write.table("clusters/cluster_1_chr23.bed", sep = "\t", quote = F, col.names = F, row.names = F)

chr23 %>%
  filter(V11 %in% data[data$cluster == "cluster_2", "gene"]) %>%
  select(V1, V2, V3, V4, V12, V13, V11) %>%
  write.table("clusters/cluster_2_chr23.bed", sep = "\t", quote = F, col.names = F, row.names = F)


chr23 %>%
  filter(V11 %in% data[data$cluster == "cluster_3", "gene"]) %>%
  select(V1, V2, V3, V4, V12, V13, V11) %>%
  write.table("clusters/cluster_3_chr23.bed", sep = "\t", quote = F, col.names = F, row.names = F)
