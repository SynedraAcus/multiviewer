#! /usr/bin/env python3

from argparse import ArgumentParser
from backend import parse_gff_line

parser = ArgumentParser()
parser.add_argument('-i', type=str, help='Gene ID list')
parser.add_argument('-f', type=str, help='Genomic assembly, FASTA')
# Assumes that all genes in ID list are present here. Other genes, if any, are
# ignored.
parser.add_argument('-g', type=str, help='Gene models, GFF')
args = parser.parse_args()

# Read the ID list

gene_ids = set()
for line in open(args.i):
    gene_ids.append(line.rstrip())

features = {x: [] for x in gene_ids}
# Read the GFF file
for record in open(args.g):
    feat = parse_gff_line()
    prefix = feat.get_id_prefix()
    if prefix in gene_ids:
        features[prefix].append(record)
