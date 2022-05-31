#!/usr/bin/env python3

import argparse
import sys
from typing import Tuple

import pandas as pd

BED_FIELDS = ["chr", "start", "end", "name", "score", "strand"]

def argument_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input")
    parser.add_argument("--peaks")
    parser.add_argument("--factors")
    parser.add_argument("-o", "--output")

    return parser

def get_deeptools_count_matrix(path: str) -> Tuple[list, pd.DataFrame]:
    with open(path) as count_file:
        header = count_file.readline().strip("\n#").split("\t")
    
    samples = [name.strip("'").removesuffix(".clean") for name in header[3:]]
    header = ["chr", "start", "end"] + samples
    counts = pd.read_csv(
        path, sep="\t", names=header, header=None, skiprows = 1
    )

    return samples, counts


def normalize_cpm(sample: pd.Series, size_factors: pd.DataFrame) -> pd.Series:
    return sample * size_factors.loc[sample.name, "factor"]


def main():
    parser = argument_parser()
    args = parser.parse_args()
    if len(sys.argv) < 4:
        parser.print_help()
        exit()

    # Read counts, peaks and factor tables
    samples, counts = get_deeptools_count_matrix(args.input)
    peaks = pd.read_csv(
        args.peaks, sep="\t",
        names=BED_FIELDS,
        usecols=range(4)
    )
    factors = pd.read_csv(
        args.factors, sep="\t",
        names=["sample", "reads", "factor"],
        index_col="sample",
    )
    
    # Add peak id to counts
    counts = counts.merge(peaks, on=["chr", "start", "end"])

    # Normalize counts to cpm
    counts[[f"{sample}_cpm" for sample in samples]] = counts[samples].apply(
        normalize_cpm, args=(factors,)
    )
    
    # Write output
    column_order = BED_FIELDS[:4] + samples + [f"{sample}_cpm" for sample in samples]
    counts[column_order].to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
