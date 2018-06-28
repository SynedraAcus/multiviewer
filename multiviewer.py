#! /usr/bin/env python3.6

from argparse import ArgumentParser
from backend import parse_gff_line
from sys import stderr

from blast_parser import parse_blast_file_to_hits
from multiplicates import is_duplicate

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
parser.add_argument('-b', type=str, help='BLAST TSV file. Assumed to be protein')
args = parser.parse_args()

# Read the ID list
gene_ids = set()
for line in open(args.i):
    gene_ids.add(line.rstrip())
print(f'Loaded {len(gene_ids)} gene IDs', file=stderr)

# Read the GFF file
features = {x: [] for x in gene_ids}
for record in open(args.g):
    feat = parse_gff_line(record)
    prefix = feat.get_id_prefix()
    if prefix in gene_ids:
        features[prefix].append(feat)

# IDs that are absent from the GFF data are reported and forgotten
skipped = [x for x in features if features[x] == []]
if len(skipped) > 0:
    print(f"No GFF data for the following IDs: {', '.join(skipped)}\n" +
          'These IDs will be ignored in further analysis.',
          file=stderr)

# Indexing the GFFs for FASTA traversal, read parsing and such
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
    processed = set()
    with open(args.f+'.genes', mode='w+') as gene_fasta:
        for record in SeqIO.parse(open(args.f), 'fasta'):
            if record.id in features_by_source:
                for feature in features_by_source[record.id]:
                    processed.add(feature.get_id_prefix())
                    segment = record[feature.start-args.flank_size:
                                     feature.end+args.flank_size]
                    if feature.strand == '-':
                        segment = segment.reverse_complement()
                    segment.id = feature.get_id_prefix()
                    segment.description = ''
                    SeqIO.write(segment, gene_fasta, 'fasta')
    missed = gene_ids.difference(processed)
    if len(missed) > 0:
        print(len(processed))
        print(f"No FASTA sources for the following IDs: {', '.join(missed)}\n" +
              'These genes are absent from genes FASTA, but may be used for ' +
              'read mapping if data for them is available.',
              file=stderr)

# Process the BLAST hit file
multiples_iterator = filter(is_duplicate, parse_blast_file_to_hits(args.b))
coordinate_sets = {x: [] for x in gene_ids}
for hit in multiples_iterator:
    coords = [hsp.query_pos for hsp in hit.hsps]
    coordinate_sets[hit.query_id.split('|')[2]].append(coords)
    missed = set()
    for gene in gene_ids:
        if coordinate_sets[gene] == []:
            missed.add(gene)
    if len(missed) > 0:
        print(len(processed))
        print(f"No BLAST data for the following IDs: {', '.join(missed)}\n" +
              'These genes are absent from BLAST file (or have only ' +
              'non-duplicated hits), but the read mapping will be displayed ' +
              'them, if available.',
              file=stderr)

# TODO: load SAM data for PacBio and Illumina reads


