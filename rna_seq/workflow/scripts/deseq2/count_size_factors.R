#!/usr/bin/env Rscript

library(DESeq2)

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
dds <- readRDS(args[1])
outfile <- args[2]

# Estimate count factors and fpm factors
# Adapted from DESEq2's fpm function
factors <- sizeFactors(dds)
k <- counts(dds)
fpm_factors <- 1e6 / (factors * exp(mean(log(colSums(k)))))

# Write output
m <- data.frame(colnames(dds), 1/factors, fpm_factors)
colnames(m) <- c("sample", "count_factors", "fpm_factors")
write.table(m, file = outfile, sep = "\t", quote = F, row.names = F)
