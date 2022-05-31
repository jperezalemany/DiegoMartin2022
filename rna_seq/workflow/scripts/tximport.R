#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(tximport))

# Log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Generate transcript-gene table
txdb <- makeTxDbFromGFF(snakemake@input[["gtf"]])
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Aggregate transcript counts to gene counts
txi <- tximport(
    snakemake@input[["salmon"]], type = "salmon", tx2gene = tx2gene
)
saveRDS(txi, file = snakemake@output[["txi"]])

# Name columns with sample names
colnames(txi$abundance) <- snakemake@params[["sample_names"]]

# Write TPMs to CSV file
write.csv(
    txi$abundance, file = snakemake@output[["tpm"]], 
    quote = F, col.names = NA
)