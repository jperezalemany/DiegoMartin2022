#!/usr/bin/env python3

""" prepare_annotation.py

This script creates an annotation bedfile containing exons, introns, UTRs 
and intergenic regions from a GTF.
If custom bed regions are provided, they will be substracted from intergenic regions
and included in the annotation.

Jaime Perez Alemany
23/04/2022
"""

import argparse
import sys

from pybedtools import BedTool
import pandas as pd


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
        usage="prepare_annotation_bed.py --gtf <gtf> --genome <chromSizes> --custom <bed> --output <bed>",
        description="This script creates an annotation bedfile containing exons, introns, UTRs "
        "and intergenic regions from a GTF. If custom bed regions are provided, "
        "they will be substracted fromintergenic regions and included in the annotation.",
    )
    parser.add_argument("--gtf", metavar="<gtf>", help="annotation in GTF format")
    parser.add_argument(
        "--genome",
        metavar="<chromSizes>",
        help="genome file with two columns: chromosome name and size, separated by TAB",
    )
    parser.add_argument(
        "--custom", metavar="<bed>", help="custom bed file with regions to include"
    )
    parser.add_argument("-t", "--gene-types", metavar="<path>")
    parser.add_argument(
        "-o", "--output", metavar="<bed>", help="path where to save the annotation file"
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


def get_genes_from_gtf(gtf: str, gene_types: str) -> pd.DataFrame:
    gtf = pd.read_csv(gtf, sep="\t", names=GTF_FIELDS)
    gtf["gene_id"] = gtf.attributes.apply(get_gene_id_from_gtf_attrs)
    gene_types = pd.read_csv(gene_types, sep="\t", names=["name", "type"])

    bed = (
        gtf[gtf.feature.isin(GENE_FEATURES)]
        .rename({"gene_id": "name"}, axis=1)
        .reindex(columns=BED_FIELDS)
        .merge(gene_types, on="name")
        .sort_values(["chrom", "start"])
    )

    return BedTool.from_dataframe(bed)


def gtf_to_bed(gtf: str, gene_types: str) -> BedTool:
    """Transforms a GTF into a bed at gene-level"""
    gtf = pd.read_csv(gtf, sep="\t", names=GTF_FIELDS)
    gene_types = pd.read_csv(gene_types, sep="\t", names=["gene_name", "gene_type"])

    bed = (
        gtf.assign(gene_name=gtf.attributes.apply(get_gene_id_from_gtf_attrs))
        .rename({"feature": "name"}, axis=1)
        .reindex(columns=BED_FIELDS + ["gene_name"])
        .drop_duplicates()
        .merge(gene_types, on="gene_name")
    )

    return BedTool.from_dataframe(bed)


def get_gene_body_annotation(feature_bed: BedTool) -> BedTool:
    """Returns non-redundant exons, introns and utrs from a bedfile"""
    # genes
    genes = feature_bed.filter(lambda x: x[3] in GENE_FEATURES)

    # UTRs
    utrs = BedTool.from_dataframe(
        feature_bed.filter(lambda x: "UTR" in x[3])
        .saveas()
        .to_dataframe()
        .drop_duplicates()
    )

    # exons
    exons = BedTool.from_dataframe(
        feature_bed.filter(
            lambda x: x[3] == "exon"
            or x[3] == "pseudogenic_exon"
            or x[3] == "miRNA_primary_transcript"
        )
        .saveas()
        .sort()
        .to_dataframe()
        .drop_duplicates()
        .assign(name="exon")
    )

    coding_exons = exons.filter(lambda x: x[7] == "protein_coding").saveas()
    coding_genes = genes.filter(lambda x: x[7] == "protein_coding").saveas()

    # introns
    introns = BedTool.from_dataframe(
        coding_genes.subtract(coding_exons.cat(utrs, postmerge=False))
        .to_dataframe()
        .drop_duplicates()
        .assign(name="intron")
    )

    return BedTool.cat(exons, introns, utrs, postmerge=False).sort()


def genome_bed(chromsizes: str) -> BedTool:
    """Returns a genome bedtool from a file of chromsome sizes"""
    genome = (
        pd.read_csv(chromsizes, sep="\t", names=["chrom", "end"])
        .assign(start=0)
        .reindex(columns=["chrom", "start", "end"])
    )

    return BedTool.from_dataframe(genome)


def intergenic_regions(
    genome: BedTool, genic_features: BedTool, genes: BedTool
) -> BedTool:
    """Computes intergenic regions by substracting genic features from a genome"""
    intergenic = BedTool.from_dataframe(
        genome.subtract(genic_features).to_dataframe().assign(name="intergenic")
    )
    intergenic = (
        BedTool.cat(
            intergenic.closest(genes, iu=True, D="a"),
            intergenic.closest(genes, id=True, D="a"),
            postmerge=False,
        )
        .to_dataframe(
            names=BED_FIELDS[0:4]
            + [f"gene_{field}" for field in BED_FIELDS + ["type"]]
            + ["distance"]
        )
        .reindex(
            columns=BED_FIELDS[0:4]
            + ["gene_score", "gene_strand", "gene_name", "gene_type"]
        )
    )
    intergenic = intergenic[intergenic.gene_name != "."].drop_duplicates()

    return BedTool.from_dataframe(intergenic)


def main():
    parser = argument_parser()
    args = parser.parse_args()
    if not args.gtf or not args.genome or not args.output:
        parser.print_help()
        exit()

    features = gtf_to_bed(args.gtf, args.gene_types)
    genes = get_genes_from_gtf(args.gtf, args.gene_types)
    gene_body = get_gene_body_annotation(features)

    custom = BedTool(args.custom)

    intergenic = intergenic_regions(
        genome=genome_bed(args.genome),
        genic_features=BedTool.cat(gene_body, genes, custom, postmerge=False),
        genes=genes,
    )

    annotation = (
        BedTool.cat(intergenic, gene_body, custom, postmerge=False)
        .sort()
        .to_dataframe(names=BED_FIELDS + ["gene_name", "gene_type"])
        .merge(
            genes.to_dataframe(
                usecols=range(6), names=[f"gene_{field}" for field in BED_FIELDS]
            ),
            on="gene_name",
        )
    )

    annotation[
        BED_FIELDS + [f"gene_{field}" for field in BED_FIELDS + ["type"]]
    ].to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
