#!/usr/bin/env python3

import argparse
import gzip
import os
import urllib.request


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download a FASTA file and replace its IDs with an alias",
        usage="download_fasta.py --link <http> --alias <csv> --output <fasta>"
    )
    parser.add_argument(
        "--link", metavar="<http>",
        help="HTTP link to the FASTA file"
    )
    parser.add_argument(
        "--output", metavar="<fasta>",
        help="path where the file will be stored"
    )
    parser.add_argument(
        "--alias", metavar="<csv>",
        help="two column file to replace IDs"
    )
    parser.add_argument(
        "--alias-sep", metavar="<chr>", default=",",
        help="alias file column separator"
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

    fasta = read_text_file(download_path)
    if fasta:
        os.remove(download_path)
    
    if args.alias:
        alias = read_two_column_file(args.alias, args.alias_sep)
        fasta = recursive_replace(fasta, alias)

    with open(args.output, "w") as output_path:
        output_path.write(fasta)


if __name__ == "__main__":
    main()
