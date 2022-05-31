#!/usr/bin/env python3

""" filter_peaks.py

This script filters regions included in a blacklist and/or greylist from a peak file.
If summits are provided, only peaks overlapping them are kept.

Jaime Perez Alemany
23/04/2022
"""

import argparse
import math

from pybedtools import BedTool
import pandas as pd

BED_FIELDS = ["chrom", "start", "end", "name", "score", "strand"]


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="filter_peaks.py [options] --inPeaks <bed> --outPeaks <bed>",
        description="This script filters regions included in a "
        "blacklist and/or greylist from a peak file. If summits are provided, "
        "only peaks overlapping them are kept.",
    )
    parser.add_argument(
        "--peaks",
        metavar="<path>",
        help="peaks in bed or narrow/broad peak format",
    )
    parser.add_argument(
        "--single-summits",
        metavar="<path>",
        nargs="+",
        help="peak summits in bed format",
    )
    parser.add_argument("--pooled-summits", metavar="<path>")
    parser.add_argument(
        "--output", metavar="<bed>", help="path where to store filtered summits"
    )
    parser.add_argument("-n", "--summits-per-peak", metavar="<int>", default="all")
    return parser


def get_pooled_summits(peaks: BedTool, pooled_summits: BedTool):
    return pooled_summits.intersect(peaks, wa=True, wb=True)


def get_single_summits(
    peaks: BedTool, pooled_summits: BedTool, single_summits: BedTool
):

    lost_peaks = peaks.intersect(pooled_summits, v=True)
    result = single_summits.intersect(lost_peaks, wa=True, wb=True)

    return result.sort()


def concat_sort(bedfiles: list) -> BedTool:
    """Concatenates a list of bedfiles into one sorted bedtool"""
    result = BedTool(bedfiles[0])

    for bed in bedfiles[1:]:
        result = result.cat(bed, postmerge=False)

    return result.sort()


def arrange_summit_columns(summits: BedTool, name: str) -> BedTool:
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

    peaks = BedTool(args.peaks)
    pooled_summits = BedTool(args.pooled_summits)

    name = peaks[0].name.removesuffix("_peak_1")

    summits = get_pooled_summits(peaks, pooled_summits)

    if args.single_summits:
        single_summits = get_single_summits(
            peaks, pooled_summits, concat_sort(args.single_summits)
        )
        summits = summits.cat(single_summits, postmerge=False).sort()

    summits = arrange_summit_columns(summits, name)

    if args.summits_per_peak != "all":
        summits = discard_summits(summits, int(args.summits_per_peak))

    summits.saveas(args.output)


if __name__ == "__main__":
    main()
