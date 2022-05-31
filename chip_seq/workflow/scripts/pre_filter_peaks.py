#!/usr/bin/env python3

""" filter_peaks.py

This script filters regions included in a blacklist and/or greylist from a peak file.

Jaime Perez Alemany
23/04/2022
"""

import argparse

from pybedtools import BedTool
import pandas as pd


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="filter_peaks.py [options] --inPeaks <bed> --outPeaks <bed>",
        description="This script filters regions included in a "
        "blacklist and/or greylist from a peak file. If summits are provided, "
        "only peaks overlapping them are kept.",
    )
    parser.add_argument(
        "--input",
        metavar="<path>",
        help="peaks in bed or narrow/broad peak format",
    )
    parser.add_argument(
        "--blacklist",
        metavar="<path>",
        help="bedfile with regions to exclude from peaks (i.e. mitochondrial chromosomes)",
    )
    parser.add_argument(
        "--greylist",
        metavar="<path>",
        help="bedfile with regions detected as enriched in controls of ChIPs",
    )
    parser.add_argument("--summits")
    parser.add_argument(
        "--output",
        metavar="<path>",
        help="path where to store filtered peaks",
    )

    return parser


def intersect(peaks: BedTool, regions: BedTool, *args, **kwargs) -> BedTool:
    """Wraps bedtools intersect"""
    return peaks.intersect(regions, *args, **kwargs)


def main():
    parser = argument_parser()
    args = parser.parse_args()
    if not args.input or not args.output:
        parser.print_help()
        exit()

    peaks = BedTool(args.input)

    filtered_peaks = peaks
    if args.blacklist:
        filtered_peaks = intersect(filtered_peaks, args.blacklist, v=True, wa=True)
    if args.greylist:
        filtered_peaks = intersect(
            filtered_peaks, args.greylist, v=True, wa=True, f=0.75
        )
    if args.summits:
        filtered_peaks = intersect(filtered_peaks, args.summits, wa=True)

    filtered_peaks.saveas(args.output)


if __name__ == "__main__":
    main()
