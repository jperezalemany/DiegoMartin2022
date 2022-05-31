#!/usr/bin/env python3

import argparse
import math

from pybedtools import BedTool
import pandas as pd

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
    parser.add_argument("--skip-chromosomes", nargs="+")
    parser.add_argument("-o", "--output")

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


def get_plus_one_nucs(table: pd.DataFrame) -> pd.DataFrame:

    plus_one = (
        table[(table.distance == 0) & (table.smt_log10pval < math.log10(0.05))]
        .sort_values("tss_distance")
        .groupby("gene_name")
        .first()
    )

    return plus_one.reset_index()


def main():
    parser = argument_parser()
    args = parser.parse_args()

    summits = BedTool(args.input)
    genes = get_genes_from_gtf(args.gtf, args.skip_chromosomes)

    annotated_summits = summits.closest(genes, D="b").to_dataframe(
        names=[f"summit_{field}" for field in BED_FIELDS]
        + NUC_FIELDS
        + [f"gene_{field}" for field in BED_FIELDS]
        + ["distance"]
    )

    # TTS - TSS distance
    fw = annotated_summits[annotated_summits.gene_strand == "+"]
    rv = annotated_summits[annotated_summits.gene_strand == "-"]

    annotated_summits["tss_distance"] = pd.NA
    annotated_summits["tts_distance"] = pd.NA

    annotated_summits.loc[fw.index,
                          "tss_distance"] = fw.summit_start - fw.gene_start
    annotated_summits.loc[fw.index,
                          "tts_distance"] = fw.summit_start - fw.gene_end

    annotated_summits.loc[rv.index,
                          "tss_distance"] = rv.gene_end - rv.summit_start
    annotated_summits.loc[rv.index,
                          "tts_distance"] = rv.gene_start - rv.summit_start

    # Plus one nucleosome

    annotated_summits = annotated_summits[
        [f"summit_{field}" for field in BED_FIELDS[:4]]
        + ["gene_score", "gene_strand"]
        + NUC_FIELDS
        + [f"gene_{field}" for field in BED_FIELDS]
        + ["tss_distance", "tts_distance"]
    ]

    annotated_summits.to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
