
from os import path

GENE_RE = "AT[0-5CM]G[0-9]+"

STAR_INDEX = [
    "chrLength.txt", "chrNameLength.txt", "chrName.txt", "chrStart.txt",
    "exonGeTrInfo.tab", "exonInfo.tab", "geneInfo.tab", 
    "Genome", "genomeParameters.txt",
    "Log.out", "SA", "SAindex", 
    "sjdbInfo.txt", "sjdbList.fromGTF.out.tab", "sjdbList.out.tab",
    "transcriptInfo.tab"
]

STAR_LOG = [
    "Log.final.out", "Log.out", "Log.progress.out"
]

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
        python3 workflow/scripts/download_reference.py \
            --link {params.url} \
            --alias {params.alias} \
            --output {output} 2> {log}
        """

ruleorder: clean_annotation > gffread_annotation > download_annotation

rule download_annotation:
    output:
        "resources/reference/{annotation}.gtf"
    params:
        url=config["reference"]["annotation_url"],
        exclude=" ".join(config["reference"]["annotation_filter"])
    log:
        "logs/download_reference/{annotation}.log"
    shell:
        """
        python3 workflow/scripts/download_reference.py \
            --link {params.url} \
            --exclude-features {params.exclude} \
            --output {output} 2> {log}
        """



rule gffread_annotation:
    input:
        "resources/reference/{annotation}.gtf"
    output:
        "resources/reference/{annotation}.gffread.gtf"
    log:
        "logs/gffread/{annotation}.gffread.log"
    conda:
        "../envs/gffread.yml"
    shell:
        "gffread -T -E {input} > {output} 2> {log}"


def get_gffread_warning_genes(wildcards):
    all_genes = []
    gene_regex = re.compile(GENE_RE)
    with open(f"logs/gffread/{wildcards.annotation}.gffread.log") as logfile:
        for line in logfile:
            genes = gene_regex.findall(line.strip())
            all_genes.extend(genes)
    
    string = "(" + "|".join(all_genes) + ")"
    return string


rule clean_annotation:
    input:
        "resources/reference/{annotation}.gffread.gtf"
    output:
        "resources/reference/{annotation}.gffread.clean.gtf"
    log:
        "logs/gffread/{annotation}.gffread.clean.log"
    params:
        regex=get_gffread_warning_genes
    shell:
        "cat {input} | grep -v -E \"{params.regex}\" > {output} 2> {log}"


rule star_index:
    input:
        fasta="resources/reference/{genome}.fasta",
        gtf="resources/reference/{annotation}.gffread.gtf"
    output:
        expand(
            "resources/star_index/{genome}_{annotation}/{star_index_files}",
            star_index_files=STAR_INDEX,
            allow_missing=True
        )
    params:
        outdir=lambda wildcards, output: path.dirname(output[0]),
        options=config["STAR"]["options"],
        read_length=config["STAR"]["read_length"] - 1
    log:
        "logs/star_index/{genome}_{annotation}.log"
    conda:
        "../envs/star.yml"
    threads: 4
    shell:
        """
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.outdir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.read_length} \
            {params.options} &> {log}
        """


rule star_align:
    input:
        read1="resources/fastq/{sample}_1.fastq.gz",
        read2="resources/fastq/{sample}_2.fastq.gz",
        index=expand(
            "resources/star_index/{genome}_{annotation}/{star_index_files}",
            genome=config["reference"]["genome_name"], 
            annotation=config["reference"]["annotation_name"],
            star_index_files=STAR_INDEX
        )
    output:
        gbam="results/bamfiles/raw/{sample}/Aligned.sortedByCoord.out.bam",
        txbam="results/bamfiles/raw/{sample}/Aligned.toTranscriptome.out.bam",
        logs=expand(
            "results/bamfiles/raw/{sample}/{star_log_files}",
            star_log_files=STAR_LOG,
            allow_missing=True
        ),
        sj="results/bamfiles/raw/{sample}/SJ.out.tab"
    log:
        "logs/star_align/{sample}.log"
    params:
        index=lambda wildcards, input: path.dirname(input.index[0]),
        prefix=lambda wildcards, output: path.dirname(output.gbam),
        options=config["STAR"]["options"]
    conda:
        "../envs/star.yml"
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.read1} {input.read2} \
            --readFilesCommand zcat \
            {params.options} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.prefix}/ &> {log}
        """


rule samtools_index:
    input:
        "results/bamfiles/raw/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "results/bamfiles/raw/{sample}/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/samtools.yml"
    log:
        "logs/samtools_index/{sample}.log"
    threads: 4
    shell:
        "samtools index {input} &> {log}"

rule samtools_view_clean:
    input:
        bam="results/bamfiles/raw/{sample}/Aligned.sortedByCoord.out.bam",
        bai="results/bamfiles/raw/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        "results/bamfiles/clean/single/{sample}.clean.bam"
    params:
        flags=config["samtools_view"]["flags"]
    log:
        "logs/samtools_view/{sample}.clean.log"
    conda:
        "../envs/samtools.yml"
    threads: 4
    shell:
        """
        samtools view -bSh \
            -@ {threads} \
            {params.flags} \
            {input.bam} > {output} 2> {log}
        """

rule samtools_index_clean:
    input:
        rules.samtools_view_clean.output
    output:
        "results/bamfiles/clean/single/{sample}.clean.bam.bai"
    log:
        "logs/samtools_index/{sample}.clean.log"
    conda:
        "../envs/samtools.yml"
    threads: 1
    shell:
        "samtools index {input} &> {log}"

rule picard_merge_clean:
    input:
        bamfiles=lambda wildcards: expand(
            "results/bamfiles/clean/single/{sample}.clean.bam",
            sample=SAMPLES[SAMPLES.condition == wildcards.condition].name
        )
    output:
        "results/bamfiles/clean/pooled/{condition}.clean.bam"
    log:
        "logs/picard_MergeSamFiles/{condition}.clean.log"
    params:
        input_flag=lambda wildcards, input: " ".join([f"-I={file}" for file in input.bamfiles])
    conda:
        "../envs/picard.yml"
    threads: 4
    shell:
        """
        picard MergeSamFiles {params.input_flag} \
            -O={output} &> {log}
        """
