from os import path

rule fastqc:
    input:
        "resources/raw_fastq/{sample}.fastq.gz"
    output:
        html="results/raw_fastq/{sample}_fastqc.html",
        zip="results/raw_fastq/{sample}_fastqc.zip"
    params:
        outdir=lambda wildcards, output: path.dirname(output.html)
    conda:
        "../envs/fastqc.yml"
    log:
        "logs/fastqc/{sample}.log"
    shell:
        "fastqc {input} -o {params.outdir} &> {log}"

rule phantompeakqualtools_single_replicates:
    input:
        "results/clean_bam/single/{chip_sample}.clean.bam"
    output:
        txt="results/qc/phantompeakqualtools/single/{chip_sample}.txt",
        pdf="results/qc/phantompeakqualtools/single/{chip_sample}.pdf",
        rdata="results/qc/phantompeakqualtools/single/{chip_sample}.RData"
    params:
        outdir=lambda wildcards, output: path.dirname(output.txt)
    log:
        "logs/phantompeakqualtools/{chip_sample}.log"
    conda:
        "../envs/phantompeakqualtools.yml"
    threads: 1
    shell:
        """
        Rscript workflow/scripts/run_spp.R \
            -c={input} \
            -odir={params.outdir} \
            -savd={output.rdata} \
            -savp={output.pdf} \
            -out={output.txt} &> {log}
        """

rule phantompeakqualtools_pooled_replicates:
    input:
        "results/clean_bam/pooled/{antibody}.clean.bam",
    output:
        txt="results/qc/phantompeakqualtools/pooled/{antibody}.txt",
        pdf="results/qc/phantompeakqualtools/pooled/{antibody}.pdf",
        rdata="results/qc/phantompeakqualtools/pooled/{antibody}.RData"
    params:
        outdir=lambda wildcards, output: path.dirname(output.txt)
    log:
        "logs/phantompeakqualtools/{antibody}.log"
    conda:
        "../envs/phantompeakqualtools.yml"
    threads: 1
    shell:
        """
        Rscript workflow/scripts/run_spp.R \
            -c={input} \
            -odir={params.outdir} \
            -savd={output.rdata} \
            -savp={output.pdf} \
            -out={output.txt} &> {log}
        """

rule samtools_stats:
    input:
        lambda wildcards: "results/bamfiles/clean/single/{condition}/{sample}.clean.bam".format(
            condition=SAMPLES[SAMPLES.name == wildcards.sample].condition.values[0],
            sample=SAMPLES[SAMPLES.name == wildcards.sample].name.values[0]
        )
    output:
        "results/qc/samtools_stats/{sample}.clean.txt"
    conda:
        "../envs/samtools.yml"
    threads: 4
    log:
        "logs/samtools_stats/{sample}.clean.log"
    shell:
        """
        samtools stats --threads {threads} {input} > {output} 2> {log}
        """