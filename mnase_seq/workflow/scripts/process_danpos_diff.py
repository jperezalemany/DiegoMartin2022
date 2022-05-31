#!/usr/bin/env python3

import argparse
import pandas as pd

BED_FIELDS = ["chrom", "start", "end", "name", "score", "strand"]
GENE_FEATURES = {"gene", "pseudogene", "transposable_element_gene"}
NUC_FIELDS = ["smt_score", "smt_log10pval", "fuzziness_score", "fuzziness_log10pval"]


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--table")
    parser.add_argument("-c", "--control")
    parser.add_argument("-t", "--treatment")
    parser.add_argument("-a", "--annotation")
    parser.add_argument("--skip-chromosomes", nargs="+")
    parser.add_argument("-o", "--output")

    return parser


def main():
    parser = argument_parser()

    args = parser.parse_args()

    # Read table
    table = pd.read_csv(args.table, sep="\t")
    table = table[~table.chr.isin(args.skip_chromosomes)]
    table["control_pos"] = table.chr + " " + table.control_smt_loca
    table["treat_pos"] = table.chr + " " + table.treat_smt_loca.astype(str)

    # Read annotation
    annotation = pd.read_csv(
        args.annotation,
        sep="\t",
        names=[f"summit_{field}" for field in BED_FIELDS]
        + NUC_FIELDS
        + [f"gene_{field}" for field in BED_FIELDS]
        + ["tss_distance", "tts_distance"],
    ).drop(
        columns=[
            f"summit_{field}" for field in ["chrom", "start", "end", "score", "strand"]
        ]
        + NUC_FIELDS
    )

    # Read nucleosome BEDs
    control_bed = pd.read_csv(
        args.control,
        sep="\t",
        usecols=[0, 1, 3],
        names=[f"control_{field}" for field in ["chrom", "start", "name"]],
    )
    control_bed["control_pos"] = (
        control_bed.control_chrom + " " + control_bed.control_start.astype(str)
    )

    treat_bed = pd.read_csv(
        args.treatment,
        sep="\t",
        usecols=[0, 1, 3],
        names=[f"treat_{field}" for field in ["chrom", "start", "name"]],
    )
    treat_bed["treat_pos"] = (
        treat_bed.treat_chrom + " " + treat_bed.treat_start.astype(str)
    )

    # Merge with bedfiles
    merged_table = (
        table.merge(control_bed, on="control_pos", how="outer")
        .merge(treat_bed, on="treat_pos", how="outer")
        .drop(
            columns=["center"]
            + [
                f"{condition}_{field}"
                for field in ["chrom", "start", "pos"]
                for condition in ["control", "treat"]
            ]
        )
    )

    # Merge with annotation
    merged_table = merged_table.merge(
        annotation, left_on="control_name", right_on="summit_name", how="outer"
    ).drop(columns="summit_name")

    # Correct shift
    merged_table.loc[merged_table.control_name.isna(), "control_smt_loca"] = pd.NA
    merged_table["control_smt_loca"] = merged_table.control_smt_loca.astype(
        pd.Int64Dtype()
    )
    fw = merged_table[merged_table.gene_strand == "+"]
    rv = merged_table[merged_table.gene_strand == "-"]

    merged_table.loc[fw.index, "treat2control_dis"] = (
        fw.treat_smt_loca - fw.control_smt_loca
    )
    merged_table.loc[rv.index, "treat2control_dis"] = (
        rv.control_smt_loca - rv.treat_smt_loca
    )
    merged_table = merged_table.sort_values(["chr", "start"])

    # Export
    danpos2_analysis = {
        "fuzziness": [column for column in table.columns if "fuzziness" in column],
        "occupancy": [
            column
            for column in table.columns
            if "smt" in column and "loca" not in column
        ],
        "shift": [
            "control_smt_loca",
            "treat_smt_loca",
            "treat2control_dis",
        ],
        "dynamism": ["diff_smt_loca"]
        + [column for column in table.columns if "point" in column],
    }

    for analysis, columns in danpos2_analysis.items():
        merged_table[["control_name", "treat_name", "gene_name"] + columns].to_csv(
            f"{args.output}.{analysis}.csv", index=False
        )


if __name__ == "__main__":
    main()
