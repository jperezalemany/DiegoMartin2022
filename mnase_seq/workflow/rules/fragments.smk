
rule sub_nucleosomal_reads:
    input:
        "results/bamfiles/clean/single/{sample}.clean.bam"
    output:
        "results/bamfiles/clean/sub_nucs/{sample}.sub_nucs.bam"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/alignmentSieve/{sample}.sub_nucs.log"
    threads: 4
    shell:
        """
        alignmentSieve -b {input} -o {output} --maxFragmentLength 100 -p {threads} &> {log}
        """