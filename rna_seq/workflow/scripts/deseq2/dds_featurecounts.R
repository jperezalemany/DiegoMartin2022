#!/usr/bin/env Rscript

library(DESeq2)

# Arguments
args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outfile <- args[2]
data <- args[3]
reference <- args[4]

# Load design
data <- read.csv(data, row.names = "name")
data$condition <- factor(data$condition)
data$condition <- relevel(data$condition, ref = reference)

# Load counts
count_files <- dir(indir, pattern = ".featureCounts$", full.names = T)
lengths <- read.table(count_files[1], header = T, sep = "\t")[, "Length"]
genes <- read.table(count_files[1], header = T, sep = "\t")[, "Geneid"]

counts <- data.frame(row.names = genes)

for (count_file in count_files) {
    counts <- cbind(
        counts, 
        sample = read.table(
            count_file, header = T, sep = "\t", row.names = "Geneid"
        )[, 6]
    )
}

colnames(counts) <- sub(".featureCounts", "", basename(count_files))
counts <- counts[rownames(data)]

# Build dataset
dds <- DESeqDataSetFromMatrix(
    countData = counts, colData = data, design = ~ condition
)
mcols(dds)$basepairs <- lengths

# Filter dataset
keep <- rowSums(counts(dds)) >= 10
print("Filtering lowly expressed genes")
print(table(keep))
dds <- dds[keep, ]

# Run analysis
dds <- DESeq(dds)
saveRDS(dds, file = outfile)

