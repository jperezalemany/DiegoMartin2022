rule nucleosome_summits:
    input:
        "results/danpos2_pooled/pooled/{condition}.Fnor.smooth.positions.ref_adjust.xls"
    output:
        "results/danpos2_processed/summits/{condition}_summits.bed"
    params:
        chroms=lambda wildcards: "(" + "|".join(config["blacklist"]["chroms"].split(" ")) + ")"
    log:
        "logs/bedfiles/{condition}_summits.log"
    shell:
        """
        cat {input} \
            | tail -n +2 \
            | grep -v -E "{params.chroms}" \
            | sort -k1,1 -k4,4n \
            | awk -v OFS="\\t" '{{ print $1, $4, $4 + 1, "{wildcards.condition}_pooled_nuc_"NR, ".", ".", $5, $6, $7, $8 }}' \
            > {output} 2> {log}
        """

rule nucleosome_beds:
    input:
        "results/danpos2_pooled/pooled/{condition}.Fnor.smooth.positions.ref_adjust.xls"
    output:
        "results/danpos2_processed/nucs/{condition}_nucs.bed"
    params:
        chroms=lambda wildcards: "(" + "|".join(config["blacklist"]["chroms"].split(" ")) + ")"
    log:
        "logs/bedfiles/{condition}_nucs.log"
    shell:
        """
        cat {input} \
            | tail -n +2 \
            | grep -v -E "{params.chroms}" \
            | sort -k1,1 -k4,4n \
            | awk -v OFS="\\t" '{{ print $1, $2, $3, "{wildcards.condition}_pooled_nuc_"NR, ".", ".", $5, $6, $7, $8 }}' \
            > {output} 2> {log}
        """


rule find_plus_one:
    input:
        nucs="results/danpos2_processed/summits/{condition}_summits.bed",
        gtf="resources/reference/{annotation}.gtf".format(
            annotation=config["reference"]["annotation_name"]
        )
    output:
        "results/danpos2_processed/plus_one/{condition}_plus_one.bed"
    params:
        chroms=config["blacklist"]["chroms"],
        distance=config["nucleosome"]["max_plus_one_distance"]
    log:
        "logs/bedfiles/{condition}_plus_one.log"
    conda:
        "../envs/pybedtools.yml"
    shell:
        """
        python3 workflow/scripts/annotate_plus_one.py \
            --input {input.nucs} \
            --gtf {input.gtf} \
            --skip-chromosomes {params.chroms} \
            --output {output} &> {log}
        """

rule nucleosome_annotation:
    input:
        bed="results/danpos2_processed/summits/{condition}_summits.bed",
        gtf="resources/reference/{annotation}.gtf".format(
            annotation=config["reference"]["annotation_name"]
        )
    output:
        "results/danpos2_processed/annotation/{condition}_summits.annotation.bed"
    log:
        "logs/bedfiles/{condition}_summist.annotation.bed"
    params:
        chroms=config["blacklist"]["chroms"],
    conda:
        "../envs/pybedtools.yml"
    threads: 4
    shell:
        """
        python3 workflow/scripts/annotate_nucs.py \
            --input {input.bed} \
            --gtf {input.gtf} \
            --skip-chromosomes {params.chroms} \
            --output {output} &> {log}
        """

rule process_diff_tables:
    input:
        table="results/danpos2_pooled/{treat}-{control}.positions.integrative.xls",
        control="results/danpos2_processed/summits/{control}_summits.bed",
        treat="results/danpos2_processed/summits/{treat}_summits.bed",
        annotation="results/danpos2_processed/annotation/{control}_summits.annotation.bed"
    output:
        tables=expand(
            "results/danpos2_processed/diff/{treat}_vs_{control}.{analysis}.csv",
            analysis=["fuzziness", "occupancy", "shift", "dynamism"],
            allow_missing=True
        )
    params:
        chroms=config["blacklist"]["chroms"],
        prefix=lambda wildcards, output: output.tables[0].removesuffix(".fuzziness.csv")
    log:
        "logs/bedfiles/{treat}_vs_{control}.log"
    conda:
        "../envs/pybedtools.yml"
    shell:
        """
        python3 workflow/scripts/process_danpos_diff.py \
            --table {input.table} \
            --control {input.control} \
            --treatment {input.treat} \
            --annotation {input.annotation} \
            --skip-chromosomes {params.chroms} \
            --output {params.prefix} &> {log}
        """

rule join_analysis:
    input:
        tables=expand(
            "results/danpos2_processed/diff/{contrast}.{analysis}.csv",
            contrast=DESIGN.treatment + "_vs_" + DESIGN.control,
            allow_missing=True
        )
    output:
        expand(
            "results/danpos2_processed/diff/{treatments}.integrated_{analysis}.csv",
            treatments="_".join(DESIGN.treatment),
            allow_missing=True
        )
    conda:
        "../envs/pybedtools.yml"
    log:
        "logs/bedfiles/integrated_{analysis}.log"
    shell:
        """
        python3 workflow/scripts/join_analysis.py \
            --input {input} \
            --output {output} &> {log}
        """

