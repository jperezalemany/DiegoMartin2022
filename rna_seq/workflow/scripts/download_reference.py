#!/usr/bin/env python3

import argparse
import gzip
import os
from typing import Callable
import urllib.request

GFF_FORMAT = [".gtf", ".gff", ".gff3"]

def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download a reference file (FASTA/GFF) from a link",
        usage="download_reference.py --link <http> --output <path>"
    )
    parser.add_argument(
        "--link", metavar="<http>",
        help="HTTP link to the file"
    )
    parser.add_argument(
        "--output", metavar="<path>",
        help="path where the file will be stored"
    )
    parser.add_argument(
        "--alias", metavar="<csv>",
        help="two column CSV file to replace chromosome IDs"
    )
    parser.add_argument(
        "--exclude-features", metavar="<feature>", nargs="+",
        help="genomic features to filter out from annotation"
    )

    return parser


def read_text_file(path: str) -> str:
    if path.endswith(".gz"):
        with gzip.open(path) as file:
            return file.read().decode("utf-8")

    with open(path) as file:
        return file.read()


def read_two_column_file(path: str, sep: str) -> dict:
    content = dict()
    with open(path) as file:
        for line in file:
            key, item = line.strip().split(sep)
            content[key] = item
    
    return content


def filter_write_column_file(
    file: str, sep: str, output: str, column: int, func: Callable
) -> str:
    with open(output, "w") as out:
        for line in file.strip().split("\n"):
            fields = line.strip().split(sep)
            
            if func(fields[column - 1]):
                out.write(line + "\n")


def recursive_replace(string: str, replacements: dict) -> str:
    for old, new in replacements.items():
        string = string.replace(old, new)
    
    return string


def main() -> None:
    parser = argument_parser()
    args = parser.parse_args()

    download_path = args.output
    if args.link.endswith(".gz"):
        download_path = f"{download_path}.gz"
    
    urllib.request.urlretrieve(args.link, download_path)

    reference = read_text_file(download_path)
    if reference:
        os.remove(download_path)
    
    if args.alias:
        alias = read_two_column_file(args.alias, sep=",")
        reference = recursive_replace(reference, alias)
    
    if any(f in download_path for f in GFF_FORMAT) and args.exclude_features:
        filter_write_column_file(
            reference, sep="\t",
            column=3, output=args.output,
            func=lambda x: x not in args.exclude_features
        )
    else:
        with open(args.output, "w") as output_path:
            output_path.write(reference)


if __name__ == "__main__":
    main()
