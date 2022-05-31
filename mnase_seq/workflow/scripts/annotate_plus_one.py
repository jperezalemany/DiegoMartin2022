#!/usr/bin/env python3

import argparse
import math

import pandas as pd
from pybedtools import BedTool


GTF_FIELDS = [
    "chrom",
    "origin",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]
BED_FIELDS = ["chrom", "start", "end", "name", "score", "strand"]
GENE_FEATURES = {"gene", "pseudogene", "transposable_element_gene"}
NUC_FIELDS = ["smt_score", "smt_log10pval",
              "fuzziness_score", "fuzziness_log10pval"]


def argument_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("--gtf")
    parser.add_argument("-o", "--output")
    parser.add_argument("--skip-chromosomes", nargs="+")

    return parser


def get_gene_id_from_gtf_attrs(attrs: str) -> str:
    """Returns the gene_id value from a string of GTF atrributes"""
    fields = attrs.split(";")

    for field in fields:
        if "gene_id" in field:
            _, raw_gene_id = field.strip().split(" ")
            gene_id = raw_gene_id.strip('"')

            return gene_id


def get_genes_from_gtf(gtf: str, exclude_chr: list) -> pd.DataFrame:
    gtf = pd.read_csv(gtf, sep="\t", names=GTF_FIELDS)
    gtf["gene_name"] = gtf.attributes.apply(get_gene_id_from_gtf_attrs)

    bed = (
        gtf[gtf.feature.isin(GENE_FEATURES) & ~gtf.chrom.isin(exclude_chr)]
        .rename({"gene_name": "name"}, axis=1)
        .reindex(columns=BED_FIELDS)
        .sort_values(["chrom", "start"])
    )

    return BedTool.from_dataframe(bed).saveas()


def get_tss(gene):

    if gene.strand == "+":
        return pd.Series([gene.start, gene.start + 1])

    if gene.strand == "-":
        return pd.Series([gene.end, gene.end + 1])


def get_tss_from_genes(genes: BedTool):

    data = genes.to_dataframe(names=BED_FIELDS + ["type"])

    data[["tss_start", "tss_end"]] = data.apply(get_tss, axis=1)

    data = data[
        ["chrom", "tss_start", "tss_end", "name",
            "score", "strand", "start", "end"]
    ].sort_values(["chrom", "tss_start"])

    return BedTool.from_dataframe(data)


def get_plus_one(tss: BedTool, nucs: BedTool) -> pd.DataFrame:

    result = (
        tss.closest(nucs, D="a", iu=True)
        .to_dataframe(
            names=[
                "chrom",
                "tss_start",
                "tss_end",
                "name",
                "score",
                "strand",
                "start",
                "end",
            ]
            + [f"nuc_{field}" for field in BED_FIELDS]
            + NUC_FIELDS
            + ["distance"]
        )
        .sort_values(["nuc_chrom", "nuc_start"])
    )
    result = result[
        (result.start <= result.nuc_start)
        & (result.end >= result.nuc_end)
        # & (result.distance <= max_distance)
    ]
    result = result[
        [f"nuc_{field}" for field in BED_FIELDS[:4]]
        + ["score", "strand"]
        + NUC_FIELDS
        + ["name", "distance"]
    ]

    return result


def main():

    parser = argument_parser()
    args = parser.parse_args()

    # Read input
    nuc_summits = BedTool(args.input)
    genes = get_genes_from_gtf(args.gtf, args.skip_chromosomes)

    # Find TSS and sig. nucleosomes
    tss = get_tss_from_genes(genes).saveas()
    sig_summits = nuc_summits.filter(lambda x: float(x[7]) < math.log10(0.05))

    # Find +1
    plus_one = get_plus_one(tss, sig_summits)

    # Filter fuzzy
    #plus_one = plus_one[plus_one.fuzziness_log10pval < math.log10(0.05)]

    # Export
    plus_one.to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
