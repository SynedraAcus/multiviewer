"""
Various procedures used in the drawing
"""


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
    # Todo: Should I split regions when in intron?
    exons = []
    for feature in features:
        if feature.feature_class == 'exon':
            exons.append((feature.start, feature.end))
    # Assuming all exons are on the same strand
    strand = feature.strand
    return [tuple((convert_single_coord(exons, x, strand) for x in y))
            for y in coord_set]
