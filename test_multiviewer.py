"""
Pytest-compatible tests for the multiviewer backend
"""

import pytest
from backend import parse_gff_line, GFF_feature


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
