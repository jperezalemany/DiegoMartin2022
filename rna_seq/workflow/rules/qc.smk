rule fastqc:
    input: 
        "results/reads/{sample}_{read}.fastq.gz"
    output: 
        html="results/qc/fastqc/{sample}_{read}_fastqc.html",
        zip="results/qc/fastqc/{sample}_{read}_fastqc.zip"
    params:
        outdir=lambda wildcards, output: path.dirname(output.html)
    conda:
        "../envs/fastqc.yml"
    log:
        "logs/fastqc/{sample}_{read}.log"
    shell: 
        "fastqc {input} -o {params.outdir} &> {log}"


rule samtool_stats_clean:
    input:
        "results/bamfiles/clean/single/{sample}.clean.bam"
    output:
        "results/qc/samtools_stats/clean/{sample}.txt"
    conda:
        "../envs/samtools.yml"
    log:
        "logs/samtools_stats/{sample}.clean.log"
    threads: 4
    shell:
        "samtools stats -@ {threads} {input} 1> {output} 2> {log}"
        
rule samtool_stats_raw:
    input:
        "results/bamfiles/raw/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "results/qc/samtools_stats/raw/{sample}.txt"
    conda:
        "../envs/samtools.yml"
    log:
        "logs/samtools_stats/{sample}.raw.log"
    threads: 4
    shell:
        "samtools stats -@ {threads} {input} 1> {output} 2> {log}"