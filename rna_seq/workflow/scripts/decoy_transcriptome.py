#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from subprocess import Popen, PIPE
from shlex import split
from snakemake.shell import shell


def decoy_file(genome_file, outdir):
    cmd1 = split(f"grep '^>' {genome_file}")
    cmd2 = split("cut -d ' ' -f1")
    cmd3 = split("sed 's/>//g'")

    with open(f"{outdir}/decoys.txt", "w") as outfile:
        p1 = Popen(cmd1, stdout=PIPE)
        p2 = Popen(cmd2, stdin=p1.stdout, stdout=PIPE)
        Popen(cmd3, stdin=p2.stdout, stdout=outfile)


def decoy_aware_transcriptome(genome_file, transcriptome_file, outdir):
    cmd = split(f"cat {transcriptome_file} {genome_file}")
    with open(f"{outdir}/gentrome.fasta", "w") as outfile:
        Popen(cmd, stdout=outfile)


def main():
    # logging
    sys.stderr = open(snakemake.log[0], "w")

    # decoy txt
    decoy_file(
        genome_file=snakemake.input.genome, 
        outdir=snakemake.params.outdir
    )

    # decoy aware fasta
    decoy_aware_transcriptome(
        genome_file=snakemake.input.genome, 
        transcriptome_file=snakemake.input.transcriptome, 
        outdir=snakemake.params.outdir
    )

if __name__ == "__main__":
    main()
