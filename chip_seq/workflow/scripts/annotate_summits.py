#!/usr/bin/env python3

import argparse
import sys

import pandas as pd
import numpy as np
from pybedtools import BedTool

BED_FIELDS = ["chrom", "start", "end", "name", "score", "strand"]


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="annotate_peaks.py -i <bed> --features <bed> --genes <bed> -o <output>"
    )
    parser.add_argument(
        "-s", "--summits", metavar="<path>", help="summits in bed format"
    )
    parser.add_argument("-p", "--peaks")
    parser.add_argument(
        "-f", "--features", metavar="<path>", help="annotations in bed format"
    )
    parser.add_argument(
        "-o", "--output", metavar="<path>", help="path where to store annotated summits"
    )

    return parser


def features_overlaping_summits(summits: BedTool, features: BedTool) -> pd.DataFrame:

    return (
        summits.intersect(features, wa=True, wb=True)
        .to_dataframe(
            names=[f"summit_{field}" for field in BED_FIELDS]
            + ["peak_name"]
            + [f"feature_{field}" for field in BED_FIELDS]
            + [f"gene_{field}" for field in BED_FIELDS]
            + ["gene_type"]
        )
        .drop(
            columns=[
                f"feature_{field}"
                for field in ["chrom", "start", "end", "strand", "score"]
            ]
        )
    )


def get_genes_from_features(features: BedTool) -> pd.DataFrame:

    return BedTool.from_dataframe(
        features.to_dataframe(
            usecols=[i for i in range(6, 13)],
            names=[f"gene_{field}" for field in BED_FIELDS + ["type"]],
        )
    )


def genes_inside_peaks(peaks: BedTool, genes: BedTool) -> pd.DataFrame:

    return peaks.intersect(genes, F=1, wa=True, wb=True).to_dataframe(
        names=[f"peak_{field}" for field in BED_FIELDS]
        + [f"gene_{field}" for field in BED_FIELDS + ["type"]]
    )


def peaks_to_summits(peaks: pd.DataFrame, summits: pd.DataFrame):

    result = (
        summits.merge(peaks, on="peak_name")
        .drop(
            columns=[
                f"peak_{field}"
                for field in ["chrom", "start", "end", "score", "strand"]
            ]
        )
        .drop_duplicates()
    )
    result["feature_name"] = "gene_inside_peak"

    return result


def tss_tts_distance(row):
    """Calculates TSS and TTS distances"""
    if row.gene_strand == "+":
        tss = row.summit_start - row.gene_start
        tts = row.summit_start - row.gene_end

    elif row.gene_strand == "-":
        tss = row.gene_end - row.summit_start
        tts = row.gene_start - row.summit_start
    else:
        return pd.Series([pd.NA, pd.NA])

    return pd.Series([tss, tts])


def main():
    # Parse arguments
    parser = argument_parser()
    args = parser.parse_args()
    if len(sys.argv) < 3:
        parser.print_help()
        exit()

    # read summits and genes
    summits = BedTool(args.summits)
    peaks = BedTool(args.peaks)
    features = BedTool(args.features)
    genes = get_genes_from_features(features)

    # Intersect summits with features
    annotated_summits = features_overlaping_summits(summits, features)

    # Short genes inside peaks
    inside_peaks = genes_inside_peaks(peaks, genes)
    inside_summits = peaks_to_summits(
        inside_peaks,
        summits.to_dataframe(
            names=[f"summit_{field}" for field in BED_FIELDS] + ["peak_name"]
        ),
    )
    inside_summits = inside_summits[
        ~(inside_summits.gene_name.isin(annotated_summits.gene_name))
    ]

    # Concatenate both
    column_order = (
        [f"summit_{field}" for field in BED_FIELDS]
        + ["peak_name"]
        + [f"gene_{field}" for field in BED_FIELDS + ["type"]]
        + ["feature_name"]
    )
    inside_summits = inside_summits[column_order]
    annotated_summits = annotated_summits[column_order]
    annotated = pd.concat([annotated_summits, inside_summits], axis=0).drop_duplicates()

    # Calculate TSS and TTS distance
    annotated[["tss_distance", "tts_distance"]] = annotated.apply(
        tss_tts_distance, axis=1
    )
    column_order = (
        [f"summit_{field}" for field in BED_FIELDS]
        + ["peak_name"]
        + [f"gene_{field}" for field in BED_FIELDS + ["type"]]
        + ["tss_distance", "tts_distance", "feature_name"]
    )

    # Write to file
    annotated[column_order].sort_values(["summit_chrom", "summit_start"]).to_csv(
        args.output, sep="\t", header=False, index=False
    )


if __name__ == "__main__":
    main()
