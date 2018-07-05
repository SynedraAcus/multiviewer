"""
Pytest-compatible tests for the multiviewer backend
"""

import pytest
from backend import parse_gff_line, GFF_feature, protein_coord_to_gene_coord, \
        convert_single_coord, reduce_coords


@pytest.fixture()
def feature_standard():
    """
    Gold standard GFF_feature
    :return:
    """
    feat = GFF_feature(feat_id='75:exon:37',
                       feature_class='exon',
                       start=83359,
                       end=84556,
                       strand='+',
                       parent='75',
                       source='scaffold00140')
    return feat


def test_feature_creation():
    feat = GFF_feature(feat_id='75:exon:37',
                       feature_class='exon',
                       start=83359,
                       end=84556,
                       strand='+',
                       parent='75',
                       source='scaffold00140')
    assert feat.id == '75:exon:37'
    assert feat.feature_class == 'exon'
    assert feat.start == 83359
    assert feat.end == 84556
    assert feat.strand == '+'
    assert feat.parent == '75'
    assert feat.source == 'scaffold00140'
    # The correct ID prefix
    assert feat.get_id_prefix() == '75'


def test_gff_parsing(feature_standard):
    test_line = 'scaffold00140\tmaker\texon\t83359\t84556\t.\t+\t.\tID=75:exon:37;Parent=75\n'
    feat = parse_gff_line(test_line)
    # Check that what produces is indeed a feature, and a correct one
    assert isinstance(feat, GFF_feature)
    assert feat.id == feature_standard.id
    assert feat.feature_class == feature_standard.feature_class
    assert feat.start == feature_standard.start
    assert feat.end == feature_standard.end
    assert feat.strand == feature_standard.strand
    assert feat.parent == feature_standard.parent
    assert feat.source == feature_standard.source


# Beware: these tests may be off by as much as 5 nucleotides
def test_single_coord_conversion():
    exons = [(1000, 3000), (9000, 10000)]
    # First exon
    assert convert_single_coord(exons, 50, strand='+') == 1148
    # Second exon
    assert convert_single_coord(exons, 950, strand='+') == 9847
    # Same on minus strand
    assert convert_single_coord(exons, 50, strand='-') == 9847
    assert convert_single_coord(exons, 950, strand='-') == 1148


def test_coordinate_conversions():
    # On a minus strand
    test_lines = [
        'scaffold00080\tmaker\tgene\t1000\t10000\t.\t-\t.\tID=23-gene;Name=23-gene',
        'scaffold00080\tmaker\tmRNA\t1000\t10000\t.\t-\t.\tID=23;Parent=23-gene;Name=23;_AED=0.12;_eAED=0.13;_QI=0|0|0|1|1|1|2|0|481;Note=protein;Ontology_term=GO:0016706,GO:0005506,GO:0055114,GO:0031418;EC=EC:1.14.11',
        'scaffold00080\tmaker\texon\t1000\t3000\t.\t-\t.\tID=23:exon:58;Parent=23',
        'scaffold00080\tmaker\texon\t9000\t10000\t.\t-\t.\tID=23:exon:57;Parent=23',
        'scaffold00080\tmaker\tCDS\t1000\t3000\t.\t-\t0\tID=23:cds;Parent=23',
        'scaffold00080\tmaker\tCDS\t9000\t10000\t.\t-\t2\tID=23:cds;Parent=23']
    gene_features = [parse_gff_line(x) for x in test_lines]
    # No domains overlapping intron but multiple domains
    assert protein_coord_to_gene_coord(gene_features,
                                       [(10, 300), (350, 800)]) == \
           [(9967, 9097), (2948, 1598)]
    # Intron-overlapping domain
    assert protein_coord_to_gene_coord(gene_features, [(10, 500)]) ==  [(9967, 2498)]
    # On the plus strand
    test_lines = [
        'scaffold00080\tmaker\tgene\t1000\t10000\t.\t+\t.\tID=23-gene;Name=23-gene',
        'scaffold00080\tmaker\tmRNA\t1000\t10000\t.\t+\t.\tID=23;Parent=23-gene;Name=23;_AED=0.12;_eAED=0.13;_QI=0|0|0|1|1|1|2|0|481;Note=protein;Ontology_term=GO:0016706,GO:0005506,GO:0055114,GO:0031418;EC=EC:1.14.11',
        'scaffold00080\tmaker\texon\t1000\t3000\t.\t+\t.\tID=23:exon:58;Parent=23',
        'scaffold00080\tmaker\texon\t9000\t10000\t.\t+\t.\tID=23:exon:57;Parent=23',
        'scaffold00080\tmaker\tCDS\t1000\t3000\t.\t+\t0\tID=23:cds;Parent=23',
        'scaffold00080\tmaker\tCDS\t9000\t10000\t.\t+\t2\tID=23:cds;Parent=23']
    gene_features = [parse_gff_line(x) for x in test_lines]
    assert protein_coord_to_gene_coord(gene_features, [(10, 300), (750, 1000)]) \
            == [(1028, 1898), (9247, 9997)]
    assert protein_coord_to_gene_coord(gene_features, [(200, 800)]) == [(1598, 9397)]


def test_overlap_simplification():
    # Correct coord recalculation
    assert reduce_coords([[(18, 30), (100, 200)], [(20, 30), (110, 210)]]) ==\
           {(19, 30): 2, (105, 205): 2}
    # Skipping the hits that can't be merged
    assert reduce_coords([[(18, 30), (100, 200), (50, 80)], [(20, 30), (110, 210)]],
                         overlap_cutoff=0.9) ==\
           {(19, 30): 2, (50, 80): 1, (105, 205): 2}
    # Overlap_cutoff, should not merge the last HSP
    assert reduce_coords([[(18, 30), (100, 200), (50, 80)], [(20, 30), (110, 210)]],
                         overlap_cutoff=1) == \
           {(19, 30): 2, (50, 80): 1, (100, 200): 1, (110, 210): 1}
