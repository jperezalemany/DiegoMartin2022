
rule download_gene_types:
    # TODO
    params:
        url=config["reference"]["annotation_gene_types"]
    output:
        "resources/reference/{annotation}.gene_types.txt"
    log:
        "logs/download_reference/{annotation}.gene_types.log"
    shell:
        """
        curl -L {params.url} | grep -v "!" | tail -n +2 > {output} 2> {log}
        """


def get_single_peaks(wildcards):
    samples = CHIPS[CHIPS.antibody == wildcards.antibody].name

    return expand(
        "results/peak_calling/single/{sample}_peaks.narrowPeak",
        sample=samples
    )

rule intersect_peaks:
    input:
        get_single_peaks
    output:
        bed="results/peak_analysis/intersection/{antibody}_peaks.bed",
        table="results/peak_analysis/intersection/{antibody}_peaks.tab"
    log:
        "logs/peak_analysis/{antibody}.peak_intersection.log"
    conda:
        "../envs/pybedtools.yml"
    threads: 4
    shell:
        """
        python3 workflow/scripts/intersect_peaks.py \
            --input {input} \
            --output {output.bed} \
            --table {output.table} \
            --region-type peak &> {log}
        """

def get_filter_input(wildcards):
    samples = CHIPS[CHIPS.antibody == wildcards.antibody]

    return {
        "peaks": "results/peak_analysis/intersection/{antibody}_peaks.bed",
        "blacklist": "resources/blacklist/{genome}.blacklist.bed".format(
            genome=config["reference"]["genome_name"]
        ),
        "greylist": "results/qc/chipseq_greylist/{control}.merged_greylist.bed".format(
            control=samples.control.values[0]
        ),
        "summits": "results/peak_calling/pooled/{antibody}_pooled_summits.bed".format(
            antibody=samples.antibody.values[0]
        )
    }


rule filter_peaks:
    input:
        unpack(get_filter_input),
    output:
        peaks="results/peak_analysis/filtering/{antibody}_peaks.filter.bed",
        summits="results/peak_analysis/filtering/{antibody}_summits.filter.bed"
    conda:
        "../envs/pybedtools.yml"
    log:
        "logs/peak_analysis/{antibody}.filtering.log"
    shell:
        """
        python3 workflow/scripts/filter_peaks.py \
            --inPeaks {input.peaks} \
            --inSummits {input.summits} \
            --blacklist {input.blacklist} \
            --greylist {input.greylist} \
            --outPeaks {output.peaks} \
            --outSummits {output.summits} &> {log}
        """

rule custom_annotation:
    input:
        regions="config/custom_regions.csv",
        gtf="resources/reference/{annotation}.gtf".format(
            annotation=config["reference"]["annotation_name"]
        ),
        genome="resources/reference/{genome}.chrom_sizes".format(
            genome=config["reference"]["genome_name"]
        ),
        types="resources/reference/{annotation}.gene_types.txt".format(
            annotation=config["reference"]["annotation_name"]
        )
    output:
        "resources/annotation/custom_regions.bed"
    log:
        "logs/peak_analysis/custom_regions.log"
    conda:
        "../envs/pybedtools.yml"
    threads: 4
    shell:
        """
        python3 workflow/scripts/custom_annotation.py \
            --gtf {input.gtf} \
            --regions {input.regions} \
            --types {input.types} \
            --genome {input.genome} \
            --output {output} &> {log}
        """


rule prepare_annotation:
    input:
        gtf="resources/reference/{annotation}.gtf",
        types="resources/reference/{annotation}.gene_types.txt",
        genome="resources/reference/{genome}.chrom_sizes".format(
            genome=config["reference"]["genome_name"]
        ),
        custom="resources/annotation/custom_regions.bed",

    output:
        "resources/annotation/{annotation}.bed"
    conda:
        "../envs/pybedtools.yml"
    log:
        "logs/peak_analysis/{annotation}.prepare_annotation.log"
    threads: 4
    shell:
        """
        python3 workflow/scripts/prepare_annotation.py \
            --gtf {input.gtf} \
            --gene-types {input.types} \
            --genome {input.genome} \
            --custom {input.custom} \
            --output {output} &> {log}
        """


rule annotate_summits:
    input:
        peaks="results/peak_analysis/filtering/{antibody}_peaks.filter.bed",
        summits="results/peak_analysis/filtering/{antibody}_summits.filter.bed",
        features="resources/annotation/{annotation}.bed".format(
            annotation=config["reference"]["annotation_name"]
        )
    output:
        "results/peak_analysis/annotation/{antibody}_summits.annotation.bed"
    log:
        "logs/peak_analysis/{antibody}.annotation.log"
    conda:
        "../envs/pybedtools.yml"
    threads: 4
    shell:
        """
        python3 workflow/scripts/annotate_summits.py \
            --summits {input.summits} \
            --peaks {input.peaks} \
            --features {input.features} \
            --output {output} &> {log}
        """

rule get_targets:
    input:
        "results/peak_analysis/annotation/{antibody}_summits.annotation.bed"
    output:
        "results/peak_analysis/targets/{antibody}.targets.bed"
    params:
        tss=config["peak_analysis"]["tss_distance"],
        tts=config["peak_analysis"]["tts_distance"]
    conda:
        "../envs/pybedtools.yml"
    log:
        "logs/pybedtools/{antibody}.targets.log"
    threads: 4
    shell:
        """
        python3 workflow/scripts/get_targets.py \
            --input {input} \
            --output {output} \
            --tss {params.tss} \
            --tts {params.tts}
        """