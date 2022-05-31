from os import path

rule bin_counts:
    input:
        bw=expand(
            "results/bigwigs/single/{sample}.Fnor.smooth.bw",
            sample=SAMPLES.name
        ),
        blacklist=expand(
            "results/qc/chipseq_greylist/{control}.merged_greylist.bed",
            control=DESIGN.control
        )
    output:
        "results/qc/multiBigwigSummary/counts.npz"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/multiBigwigSummary/counts.log"
    params:
        chroms=config["blacklist"]["chroms"],
        labels=lambda wildcards, input: [name.removesuffix(".Fnor.smooth.bw") for name in [path.basename(file) for file in input.bw]]
    shell:
        """
        multiBigwigSummary bins \
            -b {input.bw} \
            -l {params.labels} \
            --binSize 1000 \
            --chromosomesToSkip {params.chroms} \
            --blackListFileName {input.blacklist} \
            -p {threads} \
            -o {output} &> {log}
        """