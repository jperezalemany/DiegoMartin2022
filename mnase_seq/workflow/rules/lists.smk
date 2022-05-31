
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

def get_grouped_input(wildcards):

    sample = SAMPLES[SAMPLES.name == wildcards.sample]

    return expand(
        "results/bamfiles/clean/single/{condition}/{sample}.clean.bam",
        condition=sample.condition,
        allow_missing=True
    )

rule greylist:
    input:
        get_grouped_input
    output:
        bed="results/qc/chipseq_greylist/{sample}.clean-grey.bed"
    params:
        outdir=lambda wildcards, output: path.dirname(output.bed)
    log:
        "logs/chipseq_greylist/{sample}.bed"
    conda:
        "../envs/chipseq_greylist.yml"
    threads: 4
    shell:
        "chipseq-greylist {input} --outdir {params.outdir} &> {log}"


rule merge_greylist:
    input:
        expand(
            "results/qc/chipseq_greylist/{control_sample}.clean-grey.bed",
            control_sample=SAMPLES[SAMPLES.control.isna()].name
        )
    output:
        "results/qc/chipseq_greylist/{control}.merged_greylist.bed"
    log:
        "logs/bedtools_merge/{control}.merged_greylist.log"
    conda:
        "../envs/bedtools.yml"
    threads: 4
    shell:
        """
        (cat {input} \
            | sortBed \
            | mergeBed \
            | awk -v OFS="\\t" '{{ print $1, $2, $3, "greylist_"NR }}' \
            > {output}) 2> {log}
        """

rule annotate_greylist:
    input:
        bed="results/qc/chipseq_greylist/{control}.merged_greylist.bed",
        genes="../araport11/all_genes_types.bed"
    output:
        "results/qc/chipseq_greylist/{control}.merged_annotated_greylist.bed"
    conda:
        "../envs/bedtools.yml"
    threads: 4
    shell:
        "intersectBed -a {input.bed} -b {input.genes} -F 0.5 -f 0.75 -e -wa -wb > {output}"
        

rule included_regions:
    input:
        chrom_sizes="resources/reference/{genome}.chrom_sizes",
        blacklist="resources/blacklist/{genome}.blacklist.bed",
        greylist="results/qc/chipseq_greylist/{control}.merged_greylist.bed"
    output:
        "results/qc/included_regions/{genome}.{control}.bed"
    log:
        "logs/included_regions/{genome}.{control}.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        cat {input.blacklist} {input.greylist} \
            | sortBed \
            | complementBed -i - -g {input.chrom_sizes} \
            | awk -v OFS="\\t" '{{ print $1, $2, $3, "{wildcards.control}_included_"NR }}' \
            > {output} 2> {log}
        """

rule included_regions_bg:
    input:
        bed="results/qc/included_regions/{genome}.{control}.bed",
        chromsizes="resources/reference/{genome}.chrom_sizes"
    output:
        "results/qc/included_regions/{genome}.{control}.bg"
    log:
        "logs/included_regions/{genome}.{control}.bg.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools genomecov \
            -i {input.bed} -g {input.chromsizes} -bga \
            | sortBed > {output} 2> {log} \
        """

rule included_regions_wig:
    input:
        bg="results/qc/included_regions/{genome}.{control}.bg",
        chromsizes="resources/reference/{genome}.chrom_sizes"
    output:
        bw=temp("results/qc/included_regions/{genome}.{control}.bw"),
        wig="results/qc/included_regions/{genome}.{control}.wig"
    log:
        "logs/included_regions/{genome}.{control}.wig.log"
    conda:
        "../envs/kentools.yml"
    shell:
        """
        bedGraphToBigWig {input.bg} {input.chromsizes} {output.bw};
        bigWigToWig {output.bw} {output.wig}
        """
