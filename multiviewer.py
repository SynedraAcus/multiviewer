#! /usr/bin/env python3.6

from argparse import ArgumentParser
from copy import deepcopy
import svgwrite
from sys import stderr
from simplesam import Reader
import os

from backend import parse_gff_line, reduce_coords, convert_coord_dict
from blast_parser import parse_blast_file_to_hits
from multiplicates import is_duplicate, overlap

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
parser.add_argument('-e', type=str, help='TSV output for mulitplicates',
                    default='multiplicates.tsv')
parser.add_argument('--flank-size', type=int, default=1000,
                    help='Width of the flanking region to be extracted from FASTA')
parser.add_argument('-b', type=str, help='BLAST TSV file. Assumed to be protein')
parser.add_argument('-p', type=str, help='PacBio SAM file')
parser.add_argument('-d', type=str, default='schemes', help="""
        Output directory. This directory will be created, if it doesn't exist.
        Default: 'schemes'
        """)
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
regions_of_interest = {}
for gene_id in features:
    # Using only 'gene' class features
    for feature in features[gene_id]:
        if feature.feature_class == 'gene':
            if feature.source in features_by_source:
                features_by_source[feature.source].append(feature)
                regions_of_interest[feature.source].append((feature.start,
                                                            feature.end))
            else:
                features_by_source[feature.source] = [feature]
                regions_of_interest[feature.source] = [(feature.start,
                                                        feature.end)]


# Process the BLAST hit file
print('Loading BLAST hits...', file=stderr)
multiples_iterator = filter(is_duplicate, parse_blast_file_to_hits(args.b))
coordinate_sets = {}
for hit in multiples_iterator:
    coords = [hsp.query_pos for hsp in hit.hsps]
    name = hit.query_id.split('|')[2]
    if name in coordinate_sets:
        coordinate_sets[name].append(coords)
    else:
        coordinate_sets[name] = [coords]

# Averaging BLAST hits and calculating the missing genes set
blast_missed = set()
for gene_id in gene_ids:
    if gene_id not in coordinate_sets or coordinate_sets[gene_id] == []:
        blast_missed.add(gene_id)
    else:
        coordinate_sets[gene_id] = reduce_coords(coordinate_sets[gene_id])
if len(blast_missed) > 0:
    print(f"No BLAST data for the following IDs: {', '.join(blast_missed)}\n" +
          'These genes are absent from BLAST file (or have only ' +
          'non-duplicated hits), but the read mapping will be displayed ' +
          'for them, if available.',
          file=stderr)
# Converting hits to nucleotide coordinates
nucleotide_blast = {x: [] for x in gene_ids if x not in blast_missed}
for gene_id in nucleotide_blast:
    nucleotide_blast[gene_id] =\
            convert_coord_dict(features[gene_id],
                               coordinate_sets[gene_id])

# Extract the gene sequences, if args.f is set, and export data
if args.f:
    print('Loading sequences...', file=stderr)
    # Run without Biopython, if FASTA is not supplied
    from Bio import SeqIO
    processed = set()
    with open(args.f+'.genes', mode='w+') as gene_fasta:
        tsv = open(args.e, mode='w')
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
                    try:
                        blast_regions = deepcopy(nucleotide_blast[
                                                     feature.get_id_prefix()])
                    except KeyError:
                        continue
                    top_regions = sorted(sorted(blast_regions.keys(),
                                         key=lambda x: blast_regions[x])[:2],
                                         key=lambda x: x[0])
                    SeqIO.write(segment, gene_fasta, 'fasta')
                    print('\t'.join((str(x) for x in
                                (feature.get_id_prefix(),
                                record.id,
                                # Coordinates relative to extracted seq
                                top_regions[0][0]+args.flank_size-feature.start,
                                top_regions[0][1]+args.flank_size-feature.start,
                                top_regions[1][0]+args.flank_size-feature.start,
                                top_regions[1][1]+args.flank_size-feature.start,
                                # Coordinates relative to scaffold
                                top_regions[0][0] + args.flank_size,
                                top_regions[0][1] + args.flank_size,
                                top_regions[1][0] + args.flank_size,
                                top_regions[1][1] + args.flank_size,
                                ))),
                          file=tsv)
        tsv.flush()
    # Check for missing genes
    missed = gene_ids.difference(processed)
    if len(missed) > 0:
        print(len(processed))
        print(f"No FASTA sources for the following IDs: {', '.join(missed)}\n" +
              'These genes are absent from genes FASTA, but may be used for ' +
              'read mapping if data for them is available.',
              file=stderr)

# TODO: load SAM data for Illumina reads

# Loading PacBio reads
print('Loading PacBio hits...', file=stderr)
pb_reader = Reader(open(args.p))
mapped_hits = {x: [] for x in regions_of_interest}
read_counts = {x: 0 for x in gene_ids}
for hit in pb_reader:
    if hit.rname in features_by_source:
        hit_coord = hit.coords
        for gene in features_by_source[hit.rname]:
            if overlap((gene.start, gene.end), hit_coord):
                mapped_hits[hit.rname].append(hit_coord)
                read_counts[gene.get_id_prefix()] += 1


# Producing images
if os.path.isdir(args.d):
    os.chdir(args.d)
else:
    os.mkdir(args.d)
    os.chdir(args.d)
    print(f'The directory {args.d} did not exist and was created',
          file=stderr)

for gene_id in gene_ids:
    if gene_id in blast_missed:
        continue
    # Setting coordinates
    length = 0
    exons = []
    gene = None
    for feature in features[gene_id]:
        if feature.feature_class == 'gene':
            gene = feature
        elif feature.feature_class == 'exon':
            exons.append(feature)
    height = 100 + 5*len(nucleotide_blast[gene_id]) + \
             max(nucleotide_blast[gene_id].values()) + \
             7 * read_counts[gene_id]
    drawing = svgwrite.Drawing(filename=f'{gene.id}.svg',
                               size=(f'{gene.end-gene.start+200}px',
                                     f'{height}px'))
    drawing.add(drawing.rect(insert=(0, 0), size=(gene.end - gene.start + 200,
                                                  height),
                             fill=svgwrite.rgb(0xff, 0xff, 0xff)))
    for exon in exons:
        drawing.add(drawing.rect(insert=(exon.start - gene.start + 100, 20),
                                 size=(exon.end-exon.start, 20)))
    running_height = 60
    # merged BLAST hits
    for domain in nucleotide_blast[gene_id]:
        drawing.add(drawing.rect(insert=(min(domain) - gene.start + 100,
                                         running_height),
                                 size=(abs(domain[1]-domain[0]),
                                       nucleotide_blast[gene_id][domain]),
                                 stroke=svgwrite.rgb(127, 0, 0),
                                 stroke_width=5,
                                 fill='white',
                                 fill_opacity=0))
        running_height += 7
    # PacBio reads
    running_height += 20
    for mapped_region in mapped_hits[gene.source]:
        if overlap(mapped_region, (gene.start, gene.end)):
            drawing.add(drawing.line(start=(mapped_region[0]-gene.start + 100,
                                            running_height),
                                     end=(mapped_region[1]-gene.end + 100,
                                          running_height),
                                     stroke=svgwrite.rgb(0, 0, 127),
                                     stroke_width=3
                                     ))
            running_height += 5
    drawing.save()
