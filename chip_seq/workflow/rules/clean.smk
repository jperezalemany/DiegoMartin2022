from os import path
import pandas as pd

rule samtools_view_filter:
    input:
        rules.picard_MarkDuplicates.output.bam
    output:
        "results/bamfiles/clean/single/{sample}.clean.bam"
    params:
        flags=config["samtools_view_filter"]["flags"]
    log:
        "logs/samtools_view/{sample}.clean.log"
    conda:
        "../envs/samtools.yml"
    threads: 4
    shell:
        "samtools view -@ {threads} -bh {params.flags} {input} > {output} 2> {log}"

rule samtools_index_clean:
    input:
        rules.samtools_view_filter.output
    output:
        "results/bamfiles/clean/single/{sample}.clean.bam.bai"
    log:
        "logs/samtools_index/{sample}.clean.log"
    conda:
        "../envs/samtools.yml"
    threads: 1
    shell:
        "samtools index {input} &> {log}"

rule picard_MergeSamFiles:
    input:
        lambda wildcards: expand(
            "results/bamfiles/clean/single/{chip_sample}.clean.bam",
            chip_sample = CHIPS[CHIPS.antibody == wildcards.antibody].name
        )
    output:
        "results/bamfiles/clean/pooled/{antibody}.clean.bam"
    params:
        input_flag=lambda wildcards, input: expand("I={file}", file=input)
    log:
        "logs/picard_MergeSamFiles/{antibody}.pooled.log"
    conda:
        "../envs/picard.yml"
    threads: 1
    shell:
        "picard MergeSamFiles {params.input_flag} O={output} &> {log}"

rule samtools_index_merged:
    input:
        rules.picard_MergeSamFiles.output
    output:
        "results/bamfiles/clean/pooled/{antibody}.clean.bam.bai"
    log:
        "logs/samtools_index/{antibody}.clean.log"
    conda:
        "../envs/samtools.yml"
    threads: 1
    shell:
        "samtools index {input} &> {log}"

