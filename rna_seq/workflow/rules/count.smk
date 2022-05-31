rule featureCounts:
    input:
        bam="results/bamfiles/raw/{sample}/Aligned.sortedByCoord.out.bam",
        bai="results/bamfiles/raw/{sample}/Aligned.sortedByCoord.out.bam.bai",
        gtf="resources/reference/{annotation}.gffread.gtf".format(
            annotation=config["reference"]["annotation_name"]
        )
    output:
        "results/counts/{sample}.featureCounts"
    params:
        options=config["feature_counts"]["options"]
    log:
        "logs/feature_counts/{sample}.log"
    conda:
        "../envs/feature_counts.yml"
    threads: 4
    shell:
        """
        featureCounts -p \
            -T {threads} \
            -a {input.gtf} \
            {params.options} \
            -o {output} {input.bam} &> {log}
        """
    
