"""
Various procedures used in the drawing
"""

from multiplicates import overlap, overlap_len
from collections import OrderedDict

class GFF_feature():
    """
    A GFF feature data class.
    It's only designed to work with a particular dialect of GFF and is by no
    means a general purpose GFF parser.
    """
    def __init__(self, **kwargs):
        self.id = kwargs['feat_id']
        for key in ('feature_class', 'start', 'end',
                    'strand', 'parent', 'source'):
                self.__setattr__(key, kwargs[key])

    def __str__(self):
        pass

    def get_id_prefix(self):
        """
        Return the first ':' or '-'-separated element of the ID. If neither
        char is present, return ID
        :return:
        """
        return self.id.split(':')[0].split('-')[0]


def parse_gff_line(line):
    """
    Create the GFF_feature object from a GFF line
    :param line:
    :return:
    """
    l = line.rstrip().split('\t')
    key_values = l[8].split(';')
    kv_dict = {}
    for x in key_values:
        pair = x.split('=')
        kv_dict[pair[0]] = pair[1]
    return GFF_feature(feat_id=kv_dict['ID'],
                       feature_class=l[2],
                       start=int(l[3]),
                       end=int(l[4]),
                       strand=l[6],
                       parent='Parent' in kv_dict and kv_dict['Parent'] or None,
                       source=l[0])


def convert_single_coord(exons, coord, strand='+'):
    """
    Convert a single coordinate using simplified exon representation
    Exons should be sorted by (source) coordinates, whatever their strand
    :param exons:
    :param coord:
    :return:
    """
    exon_lengths = [x[1] - x[0] + 1 for x in exons]
    # A pretty ugly, but working hack. Not sure if it's not missing a codon,
    # but whatever. I need a general feel anyway.
    if strand == '-':
        c = int(sum(exon_lengths)/3) - coord
    elif strand == '+':
        c = coord
    else:
        raise ValueError("Strand should be either \'+\' or \'-\'")
    x = 0
    for index, step in enumerate(exon_lengths):
        if x + step >= c * 3:
            return exons[index][0] - x + (c - 1) * 3 + 1
        else:
            x += step


def protein_coord_to_gene_coord(features, coord_set):
    """
    Take a set of gene features and a set of coordinates (in aminoacids),
    return these coordinates in nucleotides.

    Considers only 'exon' features for intron processing
    :param features:
    :param coord_set:
    :return:
    """
    exons = []
    for feature in features:
        if feature.feature_class == 'exon':
            exons.append((feature.start, feature.end))
    # Assuming all exons are on the same strand
    strand = feature.strand
    return [tuple((convert_single_coord(exons, x, strand) for x in y))
            for y in coord_set]


def convert_coord_dict(features, coords):
    """
    Same as `protein_coord_to_gene_coord` except that it takes an OrderedDict
    whose keys are coordinate pair and returns the OrderedDict where coordinates
    are converted and keys are preserved
    :param features:
    :param coords:
    :return:
    """
    r = OrderedDict()
    exons = []
    for feature in features:
        if feature.feature_class == 'exon':
            exons.append((feature.start, feature.end))
    strand = feature.strand
    for coord_pair in coords:
        r[(convert_single_coord(exons, coord_pair[0], strand),
           convert_single_coord(exons, coord_pair[1], strand))] =\
                    coords[coord_pair]
    return r


def reduce_coord_list(coord_set, initial_counts, overlap_cutoff):
    """
    Reduce a coordinate set, averaging overlapping coordinates.
    Coordinates are considered overlapping when at least `overlap_cutoff`
    of either pair of coordinates is overlapping with the other.
    :param coord_set:
    :param counts:
    :param overlap_cutoff:
    :return:
    """
    assert len(coord_set) == len(initial_counts)
    r = [coord_set[0]]
    counts = [initial_counts[0]]
    for hsp_index, hsp in enumerate(coord_set[1:]):
        merged = False
        for index, coord in enumerate(r):
            if overlap(coord, hsp):
                l = overlap_len(coord, hsp)
                hsp_len = hsp[1] - hsp[0]
                coord_len = coord[1] - coord[0]
                if l >= overlap_cutoff * hsp_len or \
                        l > overlap_cutoff * coord_len:
                    # Forget the smaller one if the other is at least twice
                    # bigger
                    if hsp_len >= coord_len * 2:
                        r[index] = hsp
                    elif coord_len >= hsp_len * 2:
                        break
                    # Average position if they are ~equal
                    else:
                        r[index] = round((hsp[0] + coord[0]) / 2), \
                                   round((hsp[1] + coord[1]) / 2)
                    counts[index] += initial_counts[hsp_index + 1]
                    merged = True
                    break
        if not merged:
            r.append(hsp)
            counts.append(initial_counts[hsp_index + 1])
    return OrderedDict((r[x], counts[x]) for x in range(len(r)))


def reduce_coords(coord_set, overlap_cutoff=0.8):
    """
    Wrapper around reduce_coord_list that converts BLAST hit data and runs
    a second reduction round
    :param coord_set: A list of lists of coord tuples
    :param overlap_cutoff:
    :return:
    """
    coords = [hsp for hit in coord_set for hsp in hit]
    c = [1 for x in coords]
    r = reduce_coord_list(coords, c, overlap_cutoff)
    # A se
    return reduce_coord_list(tuple(r.keys()), tuple(r.values()), overlap_cutoff)

    # r = list(coord_set[0])
    # counts = [1 for _ in r]
    # # A hit is a collection of HSP coordinate pairs
    # # Non-blast coordinates are treated similarly
    # for hit in coord_set[1:]:
    #     for hsp in hit:
    #         merged = False
    #         for index, coord in enumerate(r):
    #             if overlap(coord, hsp):
    #                 l = overlap_len(coord, hsp)
    #                 hsp_len = hsp[1] - hsp[0]
    #                 coord_len = coord[1] - coord[0]
    #                 if l >= overlap_cutoff * hsp_len or \
    #                         l > overlap_cutoff * coord_len:
    #                     # Forget the smaller one if the other is at least twice
    #                     # bigger
    #                     if hsp_len >= coord_len * 2:
    #                         r[index] = hsp
    #                     elif coord_len >= hsp_len * 2:
    #                         break
    #                     # Average position if they are ~equal
    #                     else:
    #                         r[index] = round((hsp[0] + coord[0]) / 2), \
    #                                    round((hsp[1] + coord[1]) / 2)
    #                     counts[index] += 1
    #                     merged = True
    #                     break
    #         if not merged:
    #             r.append(hsp)
    #             counts.append(1)
    # return {r[x]: counts[x] for x in range(len(r))}
