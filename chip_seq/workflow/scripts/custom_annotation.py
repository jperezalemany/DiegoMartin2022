#!/usr/bin/env python3

""" custom_annotation.py

This script creates custom annotations in bed format from a GTF file, a chromsome sizes file
and a CSV file where custom regions are defined. The latter must have seven columns:
    1. region. A name for the new region
    2. refpoint (TSS or TTS). Whether to define the region in reference to the
       TSS or the TTS of genes
    3. upstream. Number of bases up from the refpoint
    4. downstream. Number of bases down from the refpoint
    5. priority. The regions with the highest priority (1) are left intact.
       These regions are substracted from those of lower priority
    6. gene_types. GTF features of genes used to build the regions, i.e. CDS,mRNA
    7. gene_overlap. Wether overlaps with genic regions should be substracted

Example custom regions file:
region,refpoint,upstream,downstream,priority,gene_types,gene_overlap
proximal_promoter,TSS,500,0,1,CDS,False
terminator,TTS,0,500,1,CDS,False
distal_promoter,TSS,2000,-500,2,CDS,False
"""

import argparse
import sys

import pandas as pd
from pybedtools import BedTool


GTF_FIELDS = [
    "chrom",
    "origin",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]
BED_FIELDS = ["chrom", "start", "end", "name", "score", "strand"]
GENE_FEATURES = {"gene", "pseudogene", "transposable_element_gene"}


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage=(
            "custom_annotation.py --gtf <gtf> --regions <csv> "
            "--genome <chromSizes> --output <bed>"
        ),
        description=(
            "This script creates custom annotations in bed format from a GTF file, "
            "a chromsome sizes file and a CSV file where custom regions are defined."
        ),
    )
    parser.add_argument("--gtf", metavar="<gtf>", help="path to a GTF annotation file")
    parser.add_argument(
        "--regions",
        metavar="<csv>",
        help="coma-separated file with seven columns: region_name, "
        "refpoint (TSS or TTS), upstream, downstream, priority, gene_types "
        "and gene overlap (True or False)",
    )
    parser.add_argument(
        "--genome",
        metavar="<chromSizes>",
        help="tab-separated file with two columns: chromosome name and length",
    )
    parser.add_argument("--types", metavar="<txt>")
    parser.add_argument(
        "--output",
        metavar="<bed>",
        help="path where to store the new annotation in BED format",
    )

    return parser


def get_gene_id_from_gtf_attrs(attrs: str) -> str:
    """Returns the gene_id value from a string of GTF atrributes"""
    fields = attrs.split(";")

    for field in fields:
        if "gene_id" in field:
            _, raw_gene_id = field.strip().split(" ")
            gene_id = raw_gene_id.strip('"')

            return gene_id


def get_genes_from_gtf(gtf: str, gene_types: set = None) -> BedTool:
    """Returns a bed file of genes from a gtf"""
    gtf = pd.read_csv(gtf, sep="\t", names=GTF_FIELDS)
    gtf["gene_name"] = gtf.attributes.apply(get_gene_id_from_gtf_attrs)

    bed = (
        gtf[gtf.feature.isin(GENE_FEATURES)]
        .rename({"gene_name": "name"}, axis=1)
        .reindex(columns=BED_FIELDS)
    )
    if gene_types is not None:
        genes_in_types = gtf[gtf.feature.isin(gene_types)].gene_name.to_list()
        bed = bed[bed.name.isin(genes_in_types)]

    return BedTool.from_dataframe(bed)


def get_tss_coordinates(gene: pd.Series) -> pd.Series:
    """Returns the TTS coordinate of a gene"""
    if gene.strand == "+":
        gene.end = gene.start + 1

    elif gene.strand == "-":
        gene.start = gene.end - 1

    return gene


def get_tts_coordinates(gene: pd.Series) -> pd.Series:
    """Returns the TSS coordinate of a gene"""
    if gene.strand == "+":
        gene.start = gene.end
        gene.end = gene.end + 1

    elif gene.strand == "-":
        gene.end = gene.start
        gene.start = gene.end - 1

    return gene


def main():
    # Parse arguments
    parser = argument_parser()
    args = parser.parse_args()
    if len(sys.argv) < 3:
        parser.print_help()
        exit()

    # Read regions
    regions = pd.read_csv(args.regions, sep="\t").sort_values("priority")
    genes = get_genes_from_gtf(args.gtf)
    gene_types = pd.read_csv(args.types, sep="\t", names=["gene_name", "gene_type"])

    # Create new regions
    new_regions = []
    for priority in regions.priority.unique():
        new_regions_same_level = []
        for _, region in regions[regions.priority == priority].iterrows():
            if region.refpoint == "TSS":
                func = get_tss_coordinates
            elif region.refpoint == "TTS":
                func = get_tts_coordinates
            new_region = BedTool.from_dataframe(
                get_genes_from_gtf(args.gtf, gene_types=region.gene_types.split(","))
                .to_dataframe()
                .apply(func, axis=1)
                .rename({"name": "gene_name"}, axis=1)
                .assign(name=region.region)
                .reindex(columns=BED_FIELDS + ["gene_name"])
            ).slop(l=region.upstream, r=region.downstream - 1, s=True, g=args.genome)

            if not region.gene_overlap:
                new_region = new_region.subtract(genes)

            for past_level_region in new_regions:
                new_region = new_region.subtract(past_level_region)

            new_regions_same_level.append(new_region)

        new_regions.extend(new_regions_same_level)

    # Concatenate all regions and export
    BedTool.cat(*new_regions, postmerge=False).to_dataframe(
        names=BED_FIELDS + ["gene_name"]
    ).merge(gene_types, on="gene_name").sort_values(["chrom", "start"]).to_csv(
        args.output, sep="\t", index=False, header=False
    )


if __name__ == "__main__":
    main()
