
def get_size_factors_input(wildcards):
    if wildcards.sample in CHIPS.control.values:
        control_name=wildcards.sample
    else:
        control_name=CHIPS[CHIPS.name == wildcards.sample].control.values[0]

    return {
        "bam": "results/bamfiles/clean/single/{sample}.clean.bam",
        "bai": "results/bamfiles/clean/single/{sample}.clean.bam.bai",
        "regions": "results/qc/included_regions/{genome}.{control_name}.bed".format(
            genome=config["reference"]["genome_name"],
            control_name=control_name
        )
    }

rule size_factors:
    input:
        unpack(get_size_factors_input)
    output:
        "results/qc/size_factors/{sample}.txt"
    params:
        flags=config["samtools_view_filter"]["flags"]
    log:
        "logs/samtools_view_count/{sample}.log"
    conda:
        "../envs/samtools.yml"
    threads: 4
    shell:
        """
        samtools view -c -@ {threads} {params.flags} -L {input.regions} {input.bam} \
            | awk -v OFS="\\t" '{{ print "{wildcards.sample}", $1, 1000000 / $1 }}' \
            > {output} 2> {log}
        """

def get_chip_fragment(wildcards, input):
    with open(input.fragment) as f:
        values = f.readline().split("\t")[2].split(",")
    
    for value in values:
        if value != "0":
            return value

def get_genomecov_input(wildcards):
    if wildcards.sample in CHIPS.control.values:
        antibody = CHIPS[CHIPS.control == wildcards.sample].antibody.values[0]
    else:
        antibody = CHIPS[CHIPS.name == wildcards.sample].antibody.values[0]

    return {
        "bam": "results/bamfiles/clean/single/{sample}.clean.bam",
        "factor": "results/qc/size_factors/{sample}.txt",
        "chrom_sizes": "resources/reference/{genome}.chrom_sizes".format(
            genome=config["reference"]["genome_name"]
        ),
        "fragment": "results/qc/phantompeakqualtools/pooled/{antibody}.txt".format(
            antibody=antibody
        )
    }

rule bedtools_genomecov:
    input:
        unpack(get_genomecov_input)
    output:
        temp("results/bigwigs/single/{sample}.cpm.bg")
    params:
        factor=lambda wildcards, input: open(input.factor).read().strip().split("\t")[2],
        fragment=get_chip_fragment
    log:
        "logs/bedtools_genomecov/{sample}.cpm.log"
    conda:
        "../envs/bedtools.yml"
    threads: 4
    shell:
        """
        bedtools genomecov -bga \
            -ibam {input.bam} \
            -scale {params.factor} \
            -fs {params.fragment} \
            | sortBed > {output} 2> {log}
        """

rule wiggletools_mean_cpm:
    input:
        lambda wildcards: expand(
            "results/bigwigs/single/{chip_samples}.cpm.bg",
            chip_samples=CHIPS[CHIPS.antibody == wildcards.antibody].name
        ),
    output:
        temp("results/bigwigs/pooled/{antibody}.cpm.bg")
    log:
        "logs/wiggletools/{antibody}.cpm.log"
    conda:
        "../envs/wiggletools.yml"
    threads: 4
    shell:
        """
        (wiggletools mean {input.bws} | wiggletools write_bg {output} - ) 2> {log}
        """

rule bgzip_tabix_single:
    input:
        "results/bigwigs/single/{sample}.cpm.bg",
    output:
        gz="results/bigwigs/single/{sample}.cpm.bg.gz",
        tbi="results/bigwigs/single/{sample}.cpm.bg.gz.tbi"
    log:
        "logs/tabix/{sample}.cpm.log"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        (bgzip -c {input} > {output.gz} && tabix {output.gz} --preset bed) 2> {log}
        """

rule bedGraphToBigWig_cpm_single:
    input:
        bg="results/bigwigs/single/{sample}.cpm.bg",
        chrom_sizes="resources/reference/{genome}.chrom_sizes".format(
            genome=config["reference"]["genome_name"]
        )
    output:
        "results/bigwigs/single/{sample}.cpm.bw"
    log:
        "logs/bedGraphToBigWig/{sample}.cpm.log"
    conda:
        "../envs/bedGraphToBigWig.yml"
    threads: 4
    shell:
        "bedGraphToBigWig {input.bg} {input.chrom_sizes} {output} &> {log}"


rule bedGraphToBigWig_cpm_pooled:
    input:
        bg="results/bigwigs/pooled/{antibody}.cpm.bg",
        chrom_sizes="resources/reference/{genome}.chrom_sizes".format(
            genome=config["reference"]["genome_name"]
        )
    output:
        "results/bigwigs/pooled/{antibody}.cpm.bw"
    log:
        "logs/bedGraphToBigWig/{antibody}.cpm.log"
    conda:
        "../envs/bedGraphToBigWig.yml"
    threads: 4
    shell:
        "bedGraphToBigWig {input.bg} {input.chrom_sizes} {output} &> {log}"


rule deeptools_bigwigCompare:
    input:
        control="results/bigwigs/single/{control}.cpm.bw",
        chip="results/bigwigs/single/{chip}.cpm.bw"
    output:
        "results/bigwigs/single/{chip}_over_{control}.log2ratio.bw"
    log:
        "logs/bigwigCompare/{chip}_over_{control}.log2ratio.log"
    conda:
        "../envs/deeptools.yml"
    threads: 4
    shell:
        """
        bigwigCompare \
            -b1 {input.chip} -b2 {input.control} \
            --binSize 1 -p {threads} \
            -o {output} &> {log}
        """


def get_wiggletools_ratio_input(wildcards):
    chip_samples = CHIPS[CHIPS.antibody == wildcards.antibody]
    control = chip_samples.control
    return {
        "bws": expand(
            "results/bigwigs/single/{chip_sample}_over_{control}.log2ratio.bw",
            chip_sample=chip_samples.name, control=control
            ),
        "chrom_sizes": "resources/reference/{genome}.chrom_sizes".format(
            genome=config["reference"]["genome_name"]
        )
    }


rule wiggletools_mean_log2ratio:
    input:
        unpack(get_wiggletools_ratio_input)
    output:
        bw=temp("results/bigwigs/pooled/{antibody}_over_{control}.log2ratio.tmp.bw"),
        bg=temp("results/bigwigs/pooled/{antibody}_over_{control}.log2ratio.bedGraph"),
    log:
        "logs/wiggletools/{antibody}_over_{control}.log2ratio.log"
    conda:
        "../envs/wiggletools.yml"
    threads: 4
    shell:
        """
        (wiggletools mean {input.bws} \
            | wigToBigWig stdin {input.chrom_sizes} {output.bw} \
            && bigWigToBedGraph {output.bw} {output.bg}) 2> {log}
        """

rule bedGraphToBigWig_log2ratio_pooled:
    input:
        bw="results/bigwigs/pooled/{antibody}_over_{control}.log2ratio.tmp.bw",
        bg="results/bigwigs/pooled/{antibody}_over_{control}.log2ratio.bedGraph",
        chrom_sizes="resources/reference/{genome}.chrom_sizes".format(
            genome=config["reference"]["genome_name"]
        )
    output:
        "results/bigwigs/pooled/{antibody}_over_{control}.log2ratio.bw",
    log:
        "logs/wiggletools/{antibody}_over_{control}.log2ratio.log"
    conda:
        "../envs/wiggletools.yml"
    threads: 4
    shell:
        "bedGraphToBigWig {input.bg} {input.chrom_sizes} {output} &> {log}"
