
library(DESeq2)
library(dplyr)
library(ComplexHeatmap)
setwd("~/../Dropbox/Bioinformatics/RNA-seq tpf chr/")

dds <- readRDS("deseq2/dds.RDS")
deg <- read.csv("deseq2/chr_tpf.integrative_dea.counts.csv")

swi <- read.csv("~/../OneDrive/Documentos/Libro1.csv", sep = ";", header = F)
swi[1, 1] <- "AT2G46020"
swi[21, 1] <- "AT5G55040"

table <- merge(deg, swi, by.x = "TAIR", by.y = "V1") %>%
  dplyr::select(V2, tpf_log2FoldChange, chr_log2FoldChange) %>%
  tibble::column_to_rownames(var = "V2") 

brewer.pal(9, "RdBu")

col_fun <- circlize::colorRamp2(breaks = seq(-4, 4, length.out = 9),
                                colors = rev(brewer.pal(9, "RdBu")),
                                )
mat <- table[swi$V2, ] %>% t()
pheatmap(mat)
library(dendsort)
row_dend = dendsort(hclust(dist(mat)))
col_dend = dendsort(hclust(dist(t(mat))))



Heatmap(mat, col = col_fun, cluster_rows = row_dend, cluster_columns = col_dend)

rld <- rlog(dds, blind = F)

dds2 <- readRDS("deseq2/dds_salmon.RDS")
rld2 <- rlog(dds2, blind = F)

set.seed(1)
mcols(dds)$basepairs <- mcols(dds)$basepair
assay(rld) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "TAIR") %>%
  merge(swi, by.x = "TAIR", by.y = "V1") %>%
  dplyr::select(-c(TAIR)) %>%
  tibble::column_to_rownames(var = "V2") %>%
  t() %>% scale() %>% t() %>%
  Heatmap(col = col_fun, row_split = 4)
