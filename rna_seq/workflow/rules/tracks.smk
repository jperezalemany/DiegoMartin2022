

def get_count_factor(wildcards, input):
    with open(input.factor) as factor:
        for line in factor:
            fields = line.strip().split("\t")
            if fields[0] == wildcards.sample:
                return fields[1]


rule bedtools_genomecov_counts:
    input:
        bam="results/bamfiles/clean/single/{sample}.clean.bam",
        factor="results/deseq2/size_factors.txt"
    output:
        temp("results/bigwigs/single/merged_strands/{sample}.counts.bg")
    log:
        "logs/bedtools_genomecov/{sample}.counts.log"
    params:
        factor=get_count_factor
    conda:
        "../envs/bedtools.yml"
    threads: 4
    shell:
        """
        bedtools genomecov \
            -ibam {input.bam} \
            -scale {params.factor} \
            -split -bga \
            | bedtools sort > {output} 2> {log}
        """

rule wiggletools_log2_counts:
    input:
        "results/bigwigs/single/merged_strands/{sample}.counts.bg"
    output:
        temp("results/bigwigs/single/merged_strands/{sample}.log2counts.bg")
    log:
        "logs/wiggletools/{sample}.log2cpm.log"
    conda:
        "../envs/wiggletools.yml"
    shell:
        """
        (wiggletools log 2 {input} \
            | wiggletools write_bg {output} - ) &> {log}
        """


rule bedGraphToBigWig_single:
    input:
        bg="results/bigwigs/single/merged_strands/{sample}.{ext}.bg",
        chromsizes="resources/reference/{genome}.fasta.fai".format(
            genome=config["reference"]["genome_name"]
        )
    output:
        "results/bigwigs/single/merged_strands/{sample}.{ext}.bw"
    log:
        "logs/bedGraphToBigWig/{sample}.{ext}.log"
    conda:
        "../envs/bedGraphToBigWig.yml"
    shell:
        "bedGraphToBigWig {input.bg} {input.chromsizes} {output} 2> {log}"


rule bgzip_tabix_single:
    input:
        "results/bigwigs/single/merged_strands/{sample}.{ext}.bg",
    output:
        gz="results/bigwigs/single/merged_strands/{sample}.{ext}.bg.gz",
        tbi="results/bigwigs/single/merged_strands/{sample}.{ext}.bg.gz.tbi"
    log:
        "logs/tabix/{sample}.{ext}.log"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        (bgzip -c {input} > {output.gz} && tabix {output.gz} --preset bed) 2> {log}
        """


rule wiggletools_mean:
    input:
        lambda wildcards: expand(
            "results/bigwigs/single/merged_strands/{sample}.{ext}.bg",
            sample=SAMPLES[SAMPLES.condition == wildcards.condition].name,
            allow_missing=True
        )
    output:
        temp("results/bigwigs/pooled/merged_strands/{condition}.{ext}.bg")
    log:
        "logs/wiggletools/{condition}.{ext}.log"
    conda:
        "../envs/wiggletools.yml"
    threads: 4
    shell:
        "(wiggletools mean {input} | wiggletools write_bg {output} - ) &> {log}"


rule bedGraphToBigWig_pooled:
    input:
        bg="results/bigwigs/pooled/merged_strands/{condition}.{ext}.bg",
        chromsizes="resources/reference/{genome}.fasta.fai".format(
            genome=config["reference"]["genome_name"]
        )
    output:
        "results/bigwigs/pooled/merged_strands/{condition}.{ext}.bw"
    log:
        "logs/bedGraphToBigWig/{condition}.{ext}.log"
    conda:
        "../envs/bedGraphToBigWig.yml"
    threads: 4
    shell:
        "bedGraphToBigWig {input.bg} {input.chromsizes} {output} &> {log}"


rule bgzip_tabix_pooled:
    input:
        "results/bigwigs/pooled/merged_strands/{condition}.{ext}.bg",
    output:
        gz="results/bigwigs/pooled/merged_strands/{condition}.{ext}.bg.gz",
        tbi="results/bigwigs/pooled/merged_strands/{condition}.{ext}.bg.gz.tbi"
    log:
        "logs/tabix/{condition}.{ext}.log"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        (bgzip -c {input} > {output.gz} && tabix {output.gz} --preset bed) 2> {log}
        """

rule deeptools_bigwigCompare:
    input:
        control="results/bigwigs/pooled/merged_strands/{control}.counts.bw",
        treat="results/bigwigs/pooled/merged_strands/{treatment}.counts.bw"
    output:
        "results/bigwigs/pooled/log2ratio/{treatment}_over_{control}.log2ratio.bw"
    log:
        "logs/bigwigCompare/{treatment}_over_{control}.log"
    conda:
        "../envs/deeptools.yml"
    threads: 4
    shell:
        """
        bigwigCompare -b1 {input.treat} -b2 {input.control} \
            -p {threads} \
            -o {output} \
            --binSize 1 \
            &> {log}
        """
