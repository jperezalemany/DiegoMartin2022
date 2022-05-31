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


def get_tss_from_genes(genes: BedTool):

    data = genes.to_dataframe(names=BED_FIELDS).assign(
        tss_start=".", tts_start=".")
    fw = data[data.strand == "+"]
    rv = data[data.strand == "-"]

    data.loc[fw.index, ["tss_start", "tss_end"]
             ] = pd.Series([fw.start, fw.start + 1])
    data.loc[rv.index, ["tss_start", "tss_end"]
             ] = pd.Series([rv.end, rv.end + 1])

    data = data[
        ["chrom", "tss_start", "tss_end", "name",
            "score", "strand", "start", "end"]
    ].sort_values(["chrom", "tss_start"])

    return BedTool.from_dataframe(data)


def annotate_closest_gene(summits: BedTool, genes: BedTool) -> pd.DataFrame:

    annotated_summits = summits.closest(genes, D="b").to_dataframe(
        names=[f"summit_{field}" for field in BED_FIELDS]
        + NUC_FIELDS
        + [f"gene_{field}" for field in BED_FIELDS]
        + ["distance"]
    )
    annotated_summits.summit_strand = annotated_summits.gene_strand

    return annotated_summits


def calculate_tts_tss_distance(summits: pd.DataFrame) -> pd.DataFrame:

    data = summits.copy()

    fw = data[data.gene_strand == "+"]
    rv = data[data.gene_strand == "-"]

    data["tss_distance"] = pd.NA
    data["tts_distance"] = pd.NA

    data.loc[fw.index, "tss_distance"] = fw.summit_start - fw.gene_start
    data.loc[fw.index, "tts_distance"] = fw.summit_start - fw.gene_end

    data.loc[rv.index, "tss_distance"] = rv.gene_end - rv.summit_start
    data.loc[rv.index, "tts_distance"] = rv.gene_start - rv.summit_start

    return data


def get_plus_one_nucs(table: pd.DataFrame) -> pd.DataFrame:

    data = table.copy()
    data["plus"] = "."

    plus_one = (
        table[(table.distance == 0) & (table.smt_log10pval < math.log10(0.05))]
        .sort_values("tss_distance")
        .groupby("gene_name")
        .first()
    )

    data.loc[data.summit_name.isin(plus_one.summit_name), "plus"] = "+1"

    return data


def get_minus_one(summits: pd.DataFrame) -> pd.DataFrame:
    data = summits.copy()
    data["minus"] = "."

    data2 = data[(data.smt_log10pval < math.log10(1))].reset_index()

    plus_ones = data2[data2.plus == "+1"]
    fw = plus_ones[plus_ones.summit_strand == "+"]
    rv = plus_ones[plus_ones.summit_strand == "-"]

    data2.loc[fw.index - 1, "minus"] = "-1"
    data2.loc[rv.index + 1, "minus"] = "-1"

    data.loc[data.comb.isin(data2[data2.minus == "-1"].comb), "minus"] = "-1"
    print(len(data[data.minus == "-1"]))

    return data


def annotate_next_plus(summits: pd.DataFrame, reference: set, n: int):
    data = summits.copy()

    candidates = data[~data.comb.isin(reference)].sort_values(
        ["summit_chrom", "summit_start"]
    )
    ref = data[data.comb.isin(reference)].sort_values(
        ["summit_chrom", "summit_start"])

    next_plus = BedTool.from_dataframe(ref).closest(
        BedTool.from_dataframe(candidates), iu=True, D="a", s=True
    )
    next_plus = next_plus.to_dataframe(usecols=[42], names=["comb"])

    data.loc[data.comb.isin(next_plus.comb), "plus"] = f"+{n}"

    return data


def main():
    parser = argument_parser()
    args = parser.parse_args()

    summits = BedTool(args.input)
    genes = get_genes_from_gtf(args.gtf, args.skip_chromosomes)

    annotated_summits = annotate_closest_gene(summits, genes)
    annotated_summits = calculate_tts_tss_distance(annotated_summits)

    annotated_summits = get_plus_one_nucs(annotated_summits)
    annotated_summits.to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":

    main()
