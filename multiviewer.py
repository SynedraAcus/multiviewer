#! /usr/bin/env python3.6

from argparse import ArgumentParser
from backend import parse_gff_line
from sys import stderr

parser = ArgumentParser()
parser.add_argument('-i', type=str, help='Gene ID list')
# Assumes that all genes in ID list are present here. Other elements, if any, are
# ignored.
parser.add_argument('-g', type=str, help='Gene models, GFF')
parser.add_argument('-f', type=str, help="""
            Genomic assembly, FASTA. If set, the part of the sequence encoding
            genes in list, as well as '--flank-size' nucleotides flanking them,
            is extracted into a separate FASTA file. If FASTA headers do not
            correspond to those in GFF file, raises an exception. 
            """)
parser.add_argument('--flank-size', type=int, default=1000,
                    help='Width of the flanking region to be extracted from FASTA')
args = parser.parse_args()

# Read the ID list
gene_ids = set()
for line in open(args.i):
    gene_ids.add(line.rstrip())
print(f'Loaded {len(gene_ids)} gene IDs\n', file=stderr)

# Read the GFF file
features = {x: [] for x in gene_ids}
for record in open(args.g):
    feat = parse_gff_line(record)
    prefix = feat.get_id_prefix()
    if prefix in gene_ids:
        features[prefix].append(feat)
        #TODO: check that no gene_ids are skipped here and in FASTA processing

# Indexing the GFFs
features_by_source = {}
for gene_id in features:
    # Using only 'gene' class features
    for feature in features[gene_id]:
        if feature.feature_class == 'gene':
            if feature.source in features_by_source:
                features_by_source[feature.source].append(feature)
            else:
                features_by_source[feature.source] = [feature]

# Extract the gene sequences, if args.f is set
if args.f:
    # Run without Biopython, if FASTA is not supplied
    from Bio import SeqIO
    with open(args.f+'.genes', mode='w+') as gene_fasta:
        for record in SeqIO.parse(open(args.f), 'fasta'):
            if record.id in features_by_source:
                for feature in features_by_source[record.id]:
                    segment = record[feature.start-args.flank_size:
                                     feature.end+args.flank_size]
                    if feature.strand == '-':
                        segment = segment.reverse_complement()
                    segment.id = feature.get_id_prefix()
                    segment.description = ''
                    SeqIO.write(segment, gene_fasta, 'fasta')

