
from os import path

rule samtools_faidx:
    input:
        "resources/reference/{genome}.fasta"
    output:
        "resources/reference/{genome}.fasta.fai"
    log:
        "logs/samtools_faidx/{genome}.log"
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools faidx {input} &> {log}"

rule chrom_sizes:
    input:
        "resources/reference/{genome}.fasta.fai",
    output:
        "resources/reference/{genome}.chrom_sizes"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        cat {input} \
            | awk -v OFS="\\t" '{{ print $1, 0, $2 }}' \
            | sortBed \
            | cut -f1,3 > {output}
        """

rule blacklist:
    input:
        rules.chrom_sizes.output
    output:
        "resources/blacklist/{genome}.blacklist.bed"
    params:
        chroms=config["blacklist"]["chroms"].replace(" ", "|")
    shell:
        """
        cat {input} \
            | grep -E '{params.chroms}' \
            | awk -v OFS="\\t" '{{ print $1, 0, $2, "blacklist_"NR }}' \
            > {output}
        """

rule greylist:
    input:
        "results/clean_bam/{control_sample}.clean.bam"
    output:
        bed="results/qc/chipseq_greylist/{control_sample}.clean-grey.bed"
    params:
        outdir=lambda wildcards, output: path.dirname(output.bed)
    log:
        "logs/chipseq_greylist/{control_sample}.bed"
    conda:
        "../envs/chipseq_greylist.yml"
    threads: 1
    shell:
        "chipseq-greylist {input} --outdir {params.outdir} &> {log}"


rule merge_greylist:
    input:
        rules.greylist.output.bed
    output:
        "results/qc/chipseq_greylist/{control_sample}.merged_greylist.bed"
    log:
        "logs/bedtools_merge/{control_sample}.merged_greylist.log"
    conda:
        "../envs/bedtools.yml"
    threads: 1
    shell:
        """
        mergeBed -i {input} \
            | awk -v OFS="\\t" '{{ print $1, $2, $3, "greylist_"NR }}' \
            > {output} 2> {log}
        """

rule included_regions:
    input:
        chrom_sizes="resources/reference/{genome}.chrom_sizes",
        blacklist="resources/blacklist/{genome}.blacklist.bed",
        greylist="results/qc/chipseq_greylist/{control_sample}.merged_greylist.bed"
    output:
        "results/qc/included_regions/{genome}.{control_sample}.bed"
    log:
        "logs/included_regions/{genome}.{control_sample}.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        cat {input.blacklist} {input.greylist} \
            | sortBed \
            | complementBed -i - -g {input.chrom_sizes} \
            | awk -v OFS="\\t" '{{ print $1, $2, $3, "{wildcards.control_sample}_included_"NR }}' \
            > {output} 2> {log}
        """