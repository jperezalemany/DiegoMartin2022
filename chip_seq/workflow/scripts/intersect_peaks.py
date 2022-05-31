#!/usr/bin/env python3

""" intersect_bed.py

This script computes the intersection between genomic rregions.
It takes two or more bedfiles and outputs a new merged bedfile and a table with the 
overlaps of the original bedfiles.

Jaime Perez Alemany
23/04/2022
"""


import argparse
import os
import sys
import re

import pandas as pd
from pybedtools import BedTool


RM_FILE_SUFFIX = {"_peaks.narrowPeak"}
RM_NAME_SUFFIX = {"-", "_"}
MERGE_OPTIONS = {"columns": [4, 5, 6], "operations": ["collapse", "mean", "first"]}


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="This script computes the intersection between genomic regions. "
        "It takes two or more bedfiles and outputs a new merged bedfile and a table with the "
        "overlaps of the original bedfiles.",
        usage="intersect_peaks.py --input <bed> --output <bed> --table <csv>",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="<path>",
        nargs="+",
        help="paths of input bedfiles",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="<path>",
        help="path where to store the merged bedfile",
    )
    parser.add_argument(
        "-t",
        "--table",
        metavar="<csv>",
        help="path where to store the overlap table",
    )
    parser.add_argument(
        "-r",
        "--region-type",
        metavar="<str>",
        default="peak",
        help="region type (default: peak)",
    )

    return parser


def get_file_names(paths: list, rm_suffix: set = None) -> list:
    """Returns the file names from a list of paths"""
    names = [os.path.basename(file) for file in paths]

    if rm_suffix:
        for suffix in rm_suffix:
            for i, name in enumerate(names):
                names[i] = name.removesuffix(suffix)

    return names


def get_common_preffix(paths: list, rm_suffix: set = None) -> str:
    """Returns the common name from a list of paths"""
    name = os.path.commonprefix(paths)
    if rm_suffix:
        for suffix in rm_suffix:
            name = name.removesuffix(suffix)

    return name


def concat_sort(bedfiles: list) -> BedTool:
    """Concatenates a list of bedfiles into one sorted bedtool"""
    result = BedTool(bedfiles[0])

    for bed in bedfiles[1:]:
        result = result.cat(bed, postmerge=False)

    return result.sort()


def merge(bed: BedTool, operations: list, columns: list) -> BedTool:
    """Wraps bedtools merge"""
    return bed.merge(c=columns, o=operations)


def add_id(bed: BedTool, common_name: str, region_type: str) -> BedTool:
    """Adds an identifier to each region in a bedfile"""
    data = bed.to_dataframe(names=["chrom", "start", "end", "names", "score", "strand"])
    data["name"] = [f"{common_name}_{region_type}_{i}" for i in range(1, len(data) + 1)]
    data = data[["chrom", "start", "end", "name", "score", "strand", "names"]]

    return BedTool.from_dataframe(data)


def intersection_table(bed: BedTool, samples: list) -> pd.DataFrame:
    """Creates an intersection table from the output of bedtools merge"""
    data = bed.to_dataframe(usecols=[3, 6], names=["name", "names"])
    data[samples] = data.names.apply(
        lambda names: pd.Series([sample in names for sample in samples])
    )
    data["num_regions"] = data.names.apply(lambda names: len(names.split(",")))
    data["num_samples"] = data[samples].sum(axis=1)

    return data[["name", "names", "num_regions", "num_samples"] + samples]


def main():
    parser = argument_parser()
    args = parser.parse_args()

    if len(sys.argv) < 3:
        parser.print_help()
        exit()

    sample_names = get_file_names(args.input, rm_suffix=RM_FILE_SUFFIX)
    name = get_common_preffix(sample_names, rm_suffix=RM_NAME_SUFFIX)

    bedfiles = concat_sort(args.input)
    bedfiles = merge(bedfiles, **MERGE_OPTIONS)
    bedfiles = add_id(bedfiles, name, args.region_type)
    bedfiles.to_dataframe(usecols=range(0, 6)).to_csv(
        args.output, sep="\t", header=False, index=False
    )

    table = intersection_table(bedfiles, samples=sample_names)
    table.to_csv(args.table, header=True, index=False, sep="\t")


if __name__ == "__main__":
    main()
