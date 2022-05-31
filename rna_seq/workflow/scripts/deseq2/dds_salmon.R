library(tximport)
library(DESeq2)

# Arguments
args <- commandArgs(trailingOnly = TRUE)
salmon_dir <- args[1]
samples_file <- args[2]
tx_file <- args[3]
control <- args[4]
out_file <- args[5]

# Read samples csv
samples <- read.csv(samples_file)
samples$condition <- factor(samples$condition)
samples$condition <- relevel(samples$condition, ref = control)

# Read reference transcripts
txdf <- read.csv(tx_file)

# Read salmon files
files <- file.path(salmon_dir, samples$name, "quant.sf")

# Build txi
txi <- tximport(files, type="salmon", tx2gene=txdf)

# Run DESeq2
dds <- DESeqDataSetFromTximport(txi, samples, ~condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
saveRDS(dds, file = out_file)
