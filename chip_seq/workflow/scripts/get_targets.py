#!/usr/bin/env python3

import argparse
import pandas as pd


BED_FIELDS = ["chrom", "start", "end", "name", "score", "strand"]


def argument_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", metavar="<path>")
    parser.add_argument("--tss", metavar="<int>", type=int)
    parser.add_argument("--tts", metavar="<int>", type=int)
    parser.add_argument("-o", "--output", metavar="<str>")

    return parser


def filter_by_distance(peaks: pd.DataFrame, tss: int, tts: int) -> pd.DataFrame:
    genes_inside_peaks = peaks[peaks.feature == "gene_inside_peak"]
    gene_body_peaks = peaks[((peaks.tss_distance >= 0) & (peaks.tts_distance <= 0))]

    intergenic_peaks = peaks[
        ~(peaks.gene_name.isin(gene_body_peaks.gene_name))
        & ((peaks.tss_distance.abs() <= tss) | (peaks.tts_distance.abs() <= tts))
    ]

    filtered_peaks = (
        pd.concat([genes_inside_peaks, gene_body_peaks, intergenic_peaks], axis=0)
        .drop_duplicates()
        .sort_values(["summit_chrom", "summit_start"])
    )

    return filtered_peaks


def get_targets(peaks: pd.DataFrame):
    columns = [f"gene_{field}" for field in BED_FIELDS + ["type"]]

    return peaks[columns].drop_duplicates().sort_values(["gene_chrom", "gene_start"])


def filter_targets(targets: pd.DataFrame, exclude: set):
    return targets[~targets.gene_name.isin(exclude)]


def main():

    parser = argument_parser()
    args = parser.parse_args()

    # Read peaks
    peaks = pd.read_csv(
        args.input,
        sep="\t",
        names=[f"summit_{field}" for field in BED_FIELDS]
        + ["peak_name"]
        + [f"gene_{field}" for field in BED_FIELDS + ["type"]]
        + ["tss_distance", "tts_distance", "feature"],
    )

    # Filter too distant
    filtered_peaks = filter_by_distance(peaks, tss=args.tss, tts=args.tts)

    # Save all targets
    targets = get_targets(filtered_peaks)
    targets.to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
