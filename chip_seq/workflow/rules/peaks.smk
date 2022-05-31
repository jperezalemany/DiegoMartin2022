from os import path
import pandas as pd

rule effective_genome_size:
    input:
        "resources/reference/{genome}.fasta"
    output:
        "resources/reference/{genome}.effective_size.txt"
    log:
        "logs/effective_genome_size/{genome}.log"
    conda:
        "../envs/faCount.yml"
    threads: 1
    shell:
        """
        faCount {input} \
            | tail -n 1 \
            | awk '{{ print $2 - $7 }}' \
            > {output} 2> {log}
        """

def get_macs2_single_input(wildcards):
    chip = CHIPS[CHIPS.name == wildcards.chip_sample]

    return {
        "control": "results/clean_bam/single/{control}.clean.bam".format(
            control=chip.control.values[0]
        ),
        "chip": "results/clean_bam/single/{chip_sample}.clean.bam",
        "fragment": "results/qc/phantompeakqualtools/pooled/{antibody}.txt".format(
            antibody=chip.antibody.values[0]
        ),
        "genome_size": "resources/reference/{genome}.effective_size.txt".format(
            genome=config["reference"]["genome_name"]
        )
    }

def get_chip_fragment(wildcards, input):
    with open(input.fragment) as f:
        values = f.readline().split("\t")[2].split(",")
    
    for value in values:
        if value != "0":
            return value

rule macs2_callpeak_single_replicates:
    input:
        unpack(get_macs2_single_input)
    output:
        peaks="results/peak_calling/single/{chip_sample}_peaks.narrowPeak"
    params:
        outdir=lambda wildcards, output: path.dirname(output.peaks),
        genome_size=lambda wildcards, input: open(input.genome_size).read().strip(),
        fragment=get_chip_fragment
    log:
        "logs/macs2_callpeak/{chip_sample}.log"
    conda:
        "../envs/macs2.yml"
    shell:
        """
        macs2 callpeak \
            -c {input.control} \
            -t {input.chip} \
            --keep-dup all \
            -g {params.genome_size} \
            --nomodel --extsize {params.fragment} \
            -n {wildcards.chip_sample} \
            --outdir {params.outdir} &> {log}
        """

def get_macs2_pooled_input(wildcards):
    chips = CHIPS[CHIPS.antibody == wildcards.antibody]
    
    return {
        "chip": "results/clean_bam/pooled/{antibody}.clean.bam",
        "control": "results/clean_bam/single/{control}.clean.bam".format(
            control=chips.control.values[0]
        ),
        "fragment": "results/qc/phantompeakqualtools/pooled/{antibody}.txt",
        "genome_size": "resources/reference/{genome}.effective_size.txt".format(
            genome=config["reference"]["genome_name"]
        )
    }

rule macs2_callpeak_pooled_replicates:
    input:
        unpack(get_macs2_pooled_input)
    output:
        peaks="results/peak_calling/pooled/{antibody}_pooled_peaks.narrowPeak",
        summits="results/peak_calling/pooled/{antibody}_pooled_summits.bed"
    params:
        outdir=lambda wildcards, output: path.dirname(output.peaks),
        genome_size=lambda wildcards, input: open(input.genome_size).read().strip(),
        fragment=get_chip_fragment,
    log:
        "logs/macs2_callpeak/{antibody}_pooled.log"
    conda:
        "../envs/macs2.yml"
    threads: 4
    shell:
        """
        macs2 callpeak \
            -c {input.control} \
            -t {input.chip} \
            --keep-dup all \
            --call-summits \
            -g {params.genome_size} \
            --nomodel --extsize {params.fragment} \
            -n {wildcards.antibody}_pooled \
            --outdir {params.outdir} &> {log}
        """


