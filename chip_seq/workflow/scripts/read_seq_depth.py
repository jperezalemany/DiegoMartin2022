#!/usr/bin/env python3

import argparse
from pathlib import Path

import pysam

def get_read_number(bamfile, options, threads=1) -> int:
    return int(pysam.view(bamfile, "-c", options, f"-@ {threads}"))

def main():
    bamfiles = ["results/raw_bam/CHR23-7_1.markdup.bam"]
    threads = 4

    samtools_view = "-F 4 -F 1024 -q 5"
    regions = "results/qc/include_regions/TAIR10.Col0_1.bed"

    options = {
        "total_reads": "",
        "mapped_reads": "-F 4",
        "filtered_reads": samtools_view,
        "filtered_included_reads": " ".join([samtools_view, f"-L {regions}"])
    }

    header = ",".join(["name"] + list(options.keys()))

    print(header)
    for bamfile in bamfiles:
        read_numbers = [get_read_number(bamfile=bamfile, threads=threads, options=options) for options in options.values()]
        print(",".join([Path(bamfile).stem] + read_numbers))

if __name__ == "__main__":
    main()
