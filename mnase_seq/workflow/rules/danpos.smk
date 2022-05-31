
from os import path

def get_danpos_input(wildcards):

    control = SAMPLES[SAMPLES.control.isna()]
    treatments = SAMPLES[~SAMPLES.name.isin(control.name)]

    return {
        "control": expand(
            "results/bamfiles/clean/single/{path}.clean.bam",
            path=control.condition + "/" + control.name
        ),
        "treatments": expand(
            "results/bamfiles/clean/single/{path}.clean.bam",
            path=treatments.condition + "/" + treatments.name
        ),
        "regions": "results/qc/included_regions/{genome}.{control}.wig".format(
            genome=config["reference"]["genome_name"],
            control=control.condition.values[0]
        )
    }

def get_danpos_paths(wildcards, input):
    danpos_paths = []
    control_dir = path.dirname(input.control[0])
    
    treat_dirs = list(set([path.dirname(folder) for folder in input.treatments]))
    for treat_dir in treat_dirs:        
        danpos_paths.append(f"{treat_dir}/:{control_dir}/")
    
    return ",".join(danpos_paths)

rule danpos2:
    input:
        unpack(get_danpos_input)
    output:
        "results/danpos_final/reference_positions.xls"
    conda:
        "../envs/danpos2.yml"
    params:
        paths=get_danpos_paths,
    log:
        "logs/danpos2/danpos2_final.log"
    threads: 4
    shell:
        """
        python workflow/scripts/danpos2/bin/danpos.py dpos \
            {params.paths} \
            --span 1 \
            --nor_region_file {input.regions} \
            --paired 1 \
            --save 1 \
            --clonalcut 1e-10 \
            --out danpos_final &> {log}
        """