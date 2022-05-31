#!/usr/bin/env python3

import sys

from Bio import SeqIO


def merge_duplicated_transcripts(fasta: str) -> set:
    transcripts = set()
    duplicated = set()
    merged_duplicated = {}

    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in transcripts:
            duplicated.add(record.id)
        transcripts.add(record.id)

    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in duplicated:
            if record.id not in merged_duplicated:
                merged_duplicated[record.id] = record.seq
            else:
                merged_duplicated[record.id] += record.seq

    return merged_duplicated


def main():
    fasta_file = sys.argv[1]
    print(fasta_file)

    dups = merge_duplicated_transcripts(fasta_file)
    print(dups)

    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id not in dups:
            print(f">{record.id}")
            print(record.seq)

    for record, sequence in dups.items():
        print(f">{record}")
        print(sequence)


if __name__ == "__main__":
    main()
