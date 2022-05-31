#!/usr/bin/env python3

import argparse
from os import path
import re

import pandas as pd

output_regex = r"_vs_[a-z0-9]+\.[a-z]+\.csv"


def get_table_dict(paths) -> dict:
    tables = {}
    for file in paths:
        name = re.sub(output_regex, "", path.basename(file))
        tables[name] = pd.read_csv(file)

    return tables


def get_columns_to_rename(columns):
    result = []
    for column in columns:
        if column == "gene_name":
            continue
        elif column == "treat2control_dis" or "control" not in column:
            result.append(column)

    return result


def rename_column(old_name: str, new_prefix: str) -> str:
    if old_name.startswith("treat"):
        return old_name.replace("treat", new_prefix)
    else:
        return f"{new_prefix}_{old_name}"


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", nargs="+")
    parser.add_argument("-o", "--output")

    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()

    tables = get_table_dict(args.input)

    # Rename
    for name, table in tables.items():
        columns_to_rename = get_columns_to_rename(table.columns)
        mapper = {column: rename_column(column, name) for column in columns_to_rename}

        tables[name] = table.rename(columns=mapper)

    # Merge
    tables_list = list(tables.values())
    first_table = tables_list[0]
    unique_table = first_table[first_table.control_name.isna()]
    merged_table = first_table[~first_table.control_name.isna()]

    for table in tables_list[1:]:
        unique = table[table.control_name.isna()]
        common = table[~table.control_name.isna()]
        merged_table = merged_table.merge(common)
        unique_table = pd.concat([unique_table, unique])

    # Concatenate and export
    pd.concat([merged_table, unique_table]).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
