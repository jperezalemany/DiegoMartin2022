
rule filter_targets:
    input:
        "results/peak_analysis/annotation/{antibody}_summits.annotation.bed"
    output:
        "results/targets/{antibody}_targets.bed"
    conda:
        "../envs/pybedtools.yml"
    log:
        "logs/targets/{antibody}.log"
    shell:
        """
        python3 workflow/
        """
