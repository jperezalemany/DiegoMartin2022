#!/usr/bin/env python3

""" annotate_peaks.py

This script annotates peaks with their coverage over genomic features
provided as a bed file.

Jaime Perez Alemany
24/04/2022
"""

import argparse
import sys

from pybedtools import BedTool

BED_FIELDS = ["chrom", "start", "end", "name", "score", "strand"]


def argument_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(
        usage="annotate_peaks.py --input <path> --features <path> --output <path>",
        description="This script annotates peaks with their coverage over genomic features "
        "provided as a bed file.",
    )
    parser.add_argument("-i", "--input", metavar="<path>", help="input peaks")
    parser.add_argument(
        "--features", metavar="<path>", help="genomic regions to overlap with peaks"
    )
    parser.add_argument(
        "-o", "--output", metavar="<path>", help="path where to store annotated peaks"
    )

    return parser


def main():

    parser = argument_parser()
    args = parser.parse_args()

    if len(sys.argv) < 3:
        parser.print_help()
        exit()

    peaks = BedTool(args.input)

    result = peaks.intersect(args.features, wo=True).to_dataframe(
        names=BED_FIELDS
        + [f"feature_{field}" for field in BED_FIELDS]
        + [f"gene_{field}" for field in BED_FIELDS]
        + ["overlap"]
    )
    columns_to_keep = (
        BED_FIELDS
        + [f"gene_{field}" for field in BED_FIELDS]
        + ["feature_name", "overlap"]
    )
    result[columns_to_keep].to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
