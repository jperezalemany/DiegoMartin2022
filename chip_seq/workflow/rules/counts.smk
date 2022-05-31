
import numpy as np
import pandas as pd

def get_multibam_input(wildcards):
    samples = CHIPS[CHIPS.antibody == wildcards.antibody]

    return {
        "bam": expand(
            "results/bamfiles/clean/single/{sample}.clean.bam",
            sample=pd.concat([samples.name, samples.control], axis=0).drop_duplicates()
        ),
        "fragment": "results/qc/phantompeakqualtools/pooled/{antibody}.txt",
        "peaks": "results/peak_analysis/filtering/{antibody}_peaks.filter.bed"
    }

def get_chip_fragment(wildcards, input):
    with open(input.fragment) as f:
        values = f.readline().split("\t")[2].split(",")
    
    for value in values:
        if value != "0":
            return value


rule multiBamSummary:
    input:
        unpack(get_multibam_input)
    output:
        matrix="results/counts/{antibody}.npz",
        counts=temp("results/counts/{antibody}.tab")
    params:
        fragment=get_chip_fragment
    log:
        "logs/multiBamSummary/{antibody}.log"
    conda:
        "../envs/deeptools.yml"
    threads: 4
    shell:
        """
        multiBamSummary BED-file \
            -b {input.bam} \
            -o {output.matrix} \
            --BED {input.peaks} \
            --smartLabels \
            --outRawCounts {output.counts} \
            --extendReads {params.fragment} \
            -p {threads} &> {log}
        """

rule size_factors_matrix:
    input:
        expand(
            "results/qc/size_factors/{sample}.txt",
            sample=SAMPLES.name
        )
    output:
        "results/qc/size_factors/size_factors.txt"
    shell:
        "cat {input} > {output}"

rule count_matrix:
    input:
        matrix="results/counts/{antibody}.tab",
        factors="results/qc/size_factors/size_factors.txt",
        peaks="results/peak_analysis/filtering/{antibody}_peaks.filter.bed"
    output:
        "results/counts/{antibody}.csv"
    log:
        "logs/counts/{antibody}.count_matrix.log"
    conda:
        "../envs/pybedtools.yml"
    threads: 4
    shell:
        """
        python3 workflow/scripts/count_matrix.py \
            --input {input.matrix} \
            --peaks {input.peaks} \
            --factors {input.factors} \
            --output {output} &> {log}
        """

