#!/usr/bin/env python3

""" filter_peaks.py

This script filters regions included in a blacklist and/or greylist from a peak file.
If summits are provided, only peaks overlapping them are kept.

Jaime Perez Alemany
23/04/2022
"""

import argparse
import re
import sys

from pybedtools import BedTool
import pandas as pd

BED_FIELDS = ["chrom", "start", "end", "name", "score", "strand"]
BLACKLIST_INTERSECT = {"f": 1e-9, "F": 1e-9}
GREYLIST_OPTIONS = {"f": 0.75, "F": 1e-9}


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="filter_peaks.py [options] --inPeaks <bed> --outPeaks <bed>",
        description="This script filters regions included in a "
        "blacklist and/or greylist from a peak file. If summits are provided, "
        "only peaks overlapping them are kept.",
    )
    parser.add_argument(
        "--inPeaks",
        metavar="<path>",
        help="peaks in bed or narrow/broad peak format",
    )
    parser.add_argument(
        "--inSummits", metavar="<path>", help="peak summits in bed format"
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
    parser.add_argument(
        "--outPeaks",
        metavar="<path>",
        help="path where to store filtered peaks",
    )
    parser.add_argument(
        "--outSummits", metavar="<bed>", help="path where to store filtered summits"
    )
    return parser


def intersect(peaks: BedTool, regions: BedTool, *args, **kwargs) -> BedTool:
    """Wraps bedtools intersect"""
    return peaks.intersect(regions, *args, **kwargs)


def drop_duplicates(peaks: BedTool) -> BedTool:
    """Drops duplicate rows from a bedtool"""
    data = peaks.to_dataframe().drop_duplicates()
    return BedTool.from_dataframe(data)


def filter_summit_columns(summits: BedTool, name: str) -> BedTool:
    """Adds an ID to summits and returns its BED columns plus their peak ID"""
    data = summits.to_dataframe(
        names=BED_FIELDS[0:5] + [f"peak_{field}" for field in BED_FIELDS]
    )

    data = data[BED_FIELDS[0:5] + ["peak_name"]]
    data["name"] = [f"{name}_summit_{i}" for i in range(1, len(data) + 1)]
    data["strand"] = "."
    data = data[BED_FIELDS + ["peak_name"]]

    return BedTool.from_dataframe(data)


def discard_summits(summits: BedTool, summits_per_peak: int) -> BedTool:
    """Keeps n summits per peak with highest score"""
    data = summits.to_dataframe(names=BED_FIELDS + ["peak_name"])

    no_dups = data.drop_duplicates(subset="peak_name", keep=False)
    dups = data[data.duplicated(subset="peak_name", keep=False)]
    dups_filtered = dups.groupby("peak_name").apply(
        lambda summits: summits.sort_values("score", ascending=False).head(
            summits_per_peak
        )
    )

    data = (
        pd.concat([no_dups, dups_filtered])
        .sort_values(["chrom", "start"])
        .reset_index(drop=True)
    )

    return BedTool.from_dataframe(data)


def main():
    parser = argument_parser()
    args = parser.parse_args()
    if not args.inPeaks or not args.outPeaks:
        parser.print_help()
        exit()

    peaks = BedTool(args.inPeaks)
    summits = BedTool(args.inSummits) if args.inSummits else None
    name = peaks[0].name.removesuffix("_peak_1")

    filtered_peaks = peaks
    if summits:
        filtered_peaks = intersect(peaks, summits, wa=True)
    if args.blacklist:
        filtered_peaks = intersect(filtered_peaks, args.blacklist, v=True, wa=True)
    if args.greylist:
        filtered_peaks = intersect(
            filtered_peaks, args.greylist, v=True, wa=True, f=0.75
        )
    filtered_peaks = drop_duplicates(filtered_peaks)
    filtered_peaks.saveas(args.outPeaks)

    if summits:
        filtered_summits = intersect(summits, filtered_peaks, wb=True)
        filtered_summits = filter_summit_columns(filtered_summits, name=name)
        filtered_summits = discard_summits(filtered_summits, summits_per_peak=2)
        filtered_summits.saveas(args.outSummits)


if __name__ == "__main__":
    main()
