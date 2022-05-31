required_packages <- list(
    "BiocParallel",
    "DESeq2",
    "dplyr",
    "glue",
    "rtracklayer",
    "stringr"
)

# Load packages
for(package in required_packages) {
    suppressPackageStartupMessages(library(package, character.only = TRUE))
}

# Write log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Build DDS object
coldata <- read.csv(
    snakemake@input[["coldata"]], row.names = "sample_name"
    ) %>% 
    select(-c(fastq1, fastq2, strandness))
coldata[, 1:ncol(coldata)] <- lapply(coldata, as.factor)

rowdata <- import(snakemake@input[["rowdata"]]) %>%
    data.frame() %>%
    filter(type == "exon") %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    select(seqnames, start, end, width, strand, gene_id) %>%
    GRanges()

counts <- snakemake@input[["counts"]]
design <- as.formula(snakemake@params[["design"]])

if(endsWith(counts, "raw_counts.csv")) {
    counts <- read.csv(counts, row.names = 1)
    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = coldata,
        rowData = rowdata,
        design = design
    )
} else if(endsWith(counts, "txi.RDS")) {
    counts <- readRDS(counts)
    dds <- DESeqDataSetFromTximport(
        txi = counts,
        colData = coldata,
        design = design
    )
}

# Filter low counts
keep <- rowSums(counts(dds) >= 10 ) >= 3
dds <- dds[keep, ]

# Run analysis
dds <- DESeq(dds)
saveRDS(dds, file = snakemake@output[["dds"]])

# Get results
outdir <- dirname(snakemake@output[["results"]])[1]
contrasts <- snakemake@output[["results"]] %>% 
    str_remove(glue("{outdir}/")) %>%
    str_remove(".csv")

for (contrast in contrasts) {
    res <- results(
        dds, name = contrast, alpha = 0.05,
        parallel = TRUE, BPPARAM = MulticoreParam(4)
    )
    summary(res)
    if(contrast %in% resultsNames(dds)) {
        res <- lfcShrink(
            dds, coef = contrast, res = res,
            parallel = TRUE, BPPARAM = MulticoreParam(4)
        )
    }
    res %>% data.frame() %>% arrange(padj) %>% 
        write.csv(
            file = glue("{outdir}/{contrast}.csv"), quote = FALSE
        )
}
