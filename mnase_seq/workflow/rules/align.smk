
from os import path

BOWTIE2_SUFFIX = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]

rule download_genome:
    output:
        "resources/reference/{genome}.fasta"
    params:
        url=config["reference"]["genome_url"],
        alias=config["reference"]["genome_alias"]
    log:
        "logs/download_genome/{genome}.log"
    threads: 1
    shell:
        """
        python3 workflow/scripts/download_fasta.py \
            --link {params.url} \
            --alias {params.alias} \
            --output {output} 2> {log}
        """

rule bowtie2_index:
    input:
        "resources/reference/{genome}.fasta"
    output:
        index=expand(
            "resources/bowtie2_index/{genome}.{suffix}",
            suffix=BOWTIE2_SUFFIX,
            allow_missing=True
        )
    params:
        prefix=lambda wildcards, output: path.join(path.dirname(output.index[0]), wildcards.genome)
    conda:
        "../envs/bowtie2.yml"
    log:
        "logs/bowtie2_index/{genome}.log"
    threads: 4
    shell:
        "bowtie2-build --threads {threads} {input} {params.prefix} &> {log}"

rule bowtie2_align:
    input:
        index=expand(
            "resources/bowtie2_index/{genome}.{suffix}",
            genome=config["reference"]["genome_name"],
            suffix=BOWTIE2_SUFFIX
        ),
        fastq1="resources/raw_fastq/{sample}_1.fastq.gz",
        fastq2="resources/raw_fastq/{sample}_2.fastq.gz"
    
    output:
        temp("results/bamfiles/raw/{sample}.sam")
    params:
        options=config["bowtie2_align"]["options"],
        prefix=lambda wildcards, input: path.join(path.dirname(input.index[0]), config["reference"]["genome_name"])
    log:
        "logs/bowtie2_align/{sample}.log"
    conda:
        "../envs/bowtie2.yml"
    threads: 4
    shell:
        """
        bowtie2 {params.options} \
            -p {threads} \
            -x {params.prefix} \
            -1 {input.fastq1} -2 {input.fastq2} \
            > {output} 2> {log}
        """

rule samtools_view_bam:
    input:
        "results/bamfiles/raw/{sample}.sam"
    output:
        temp("results/bamfiles/raw/{sample}.bam")
    log:
        "logs/samtools_view/{sample}.log"
    conda:
        "../envs/samtools.yml"
    threads: 4
    shell:
        "samtools view -@ {threads} -bSh {input} > {output} 2> {log}"

rule samtools_sort:
    input:
        "results/bamfiles/raw/{sample}.bam"
    output:
        temp("results/bamfiles/raw/{sample}.sorted.bam")
    log:
        "logs/samtools_sort/{sample}.log"
    conda:
        "../envs/samtools.yml"
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output} {input} &> {log}"

rule picard_MarkDuplicates:
    input:
        "results/bamfiles/raw/{sample}.sorted.bam"
    output:
        bam="results/bamfiles/raw/{sample}.markdup.bam",
        txt="results/bamfiles/raw/{sample}.markdup_metrics.txt"
    log:
        "logs/picard_MarkDuplicates/{sample}.log"
    conda:
        "../envs/picard.yml"
    threads: 1
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.txt} &> {log}"

rule samtools_index_raw:
    input:
        rules.picard_MarkDuplicates.output.bam
    output:
        "results/bamfiles/raw/{sample}.markdup.bam.bai"
    log:
        "logs/samtools_index/{sample}.markdup.log"
    conda:
        "../envs/samtools.yml"
    threads: 1
    shell:
        "samtools index {input} &> {output}"
