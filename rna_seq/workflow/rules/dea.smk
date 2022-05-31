from os import path

rule deseq2_dataset_featurecounts:
    input:
        expand(
            "results/counts/{sample}.featureCounts",
            sample=SAMPLES.name
        )
    output:
        "results/deseq2/dds_featurecounts.RDS"
    params:
        indir=lambda wildcards, input: path.dirname(input[0]),
        control=config["deseq2"]["control"],
        data="config/data.csv"
    log:
        "logs/deseq2/dds_featurecounts.log"
    conda:
        "../envs/deseq2.yml"
    threads: 4
    shell:
        """
        Rscript workflow/scripts/deseq2/dds_featurecounts.R \
            {params.indir} \
            {output} \
            {params.data} \
            {params.control} \
            &> {log}
        """


rule transcripts_gene_table:
    input:
        "resources/reference/{annotation}.gffread.gtf"
    output:
        "resources/reference/{annotation}.tx2gene.csv"
    shell:
        """
        echo "TXNAME,GENEID" > {output};
        cat {input} \
            | awk -v OFS="," '{{ if ($3 == "exon") print $10, $12 }}' \
            | sed 's/[";]//g' \
            | sort | uniq >> {output}
        """



rule deseq2_dataset_salmon:
    input:
        files=expand(
            "results/salmon/alignment_mode/{sample}/quant.sf",
            sample=SAMPLES.name
        ),
        txfile="resources/reference/{annotation}.tx2gene.csv".format(
            annotation=config["reference"]["annotation_name"]
        ),
    output:
        "results/deseq2/dds_salmon.RDS"
    params:
        control=config["deseq2"]["control"],
        salmon_dir=lambda wildcards, input: "/".join(input.files[0].split("/")[:-2]),
        samples="config/data.csv"
    log:
        "logs/deseq2/dds_salmon.log"
    conda:
        "../envs/deseq2.yml"
    threads: 4
    shell:
        """
        Rscript workflow/scripts/deseq2/dds_salmon.R \
            {params.salmon_dir} \
            {params.samples} \
            {input.txfile} \
            {params.control} \
            {output} &> {log}
        """

rule deseq2_size_factors_counts:
    input:
        "results/deseq2/dds_featurecounts.RDS"
    output:
        "results/deseq2/count_size_factors.txt",
    log:
        "logs/deseq2/size_factors.log"
    conda:
        "../envs/deseq2.yml"
    threads: 4
    shell:
        """
        Rscript workflow/scripts/deseq2/count_size_factors.R \
            {input} \
            {output} &> {log}
        """
