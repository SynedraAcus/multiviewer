"""
A function to test for a multiple/nonmultiple hit.
The detection of multiplicated genes is valuable both by itself (as a chapter in
my thesis) and as a stage in a phylomic analysis. An alignment of multiplicated
gene to its non-multiplicated counterpart is unreliable, as the short one will
more or less randomly align to one of the modules of the multiple. It creates a
lot of problems: the distance calculated is one of several equally valid values
(one for each module), the alignment for ML/bayesian is needlessly bloated, a
remarkable evolutionary event is ignored altogether.
Currently I assume that a multiplication of a gene is much more common than a
previously repeat-containing gene being split in pieces. Thus, the (ancestral)
shorter form is always assumed to be in a hit, and the longer in a query.
"""

from itertools import combinations


def overlap(range1, range2):
    """
    Return True if two ranges overlap
    :param range1:
    :param range2:
    :return:
    """
    if range1[0] <= range2[1] and range2[0] <= range1[1]:
        return True
    return False


def overlap_len(range1, range2):
    """
    Return the overlap length between two ranges
    :param range1:
    :param range2:
    :return:
    """
    return min(range1[1], range2[1]) - max(range1[0], range2[0])


def is_duplicate(hit, overlap_cutoff=0.5, len_cutoff=50):
    """
    Return True if this hit covers several regions of query sequence with the
    same region of hit sequence.
    Accepts two optional parameters:
    `len_cutoff` is a minimum length (on hit) an HSP must have to be considered
    meaningful. Defaults to 50.
    `overlap_cutoff` is a minimum percentage of both hsps that should be in the
    overlap for it to be considered meaningful. Defaults to 0.5
    :param hit: BlastHit
    :param overlap_cutoff: int
    :param len_cutoff: int
    :return:
    """
    if len(hit.hsps) < 2:
        #  Obviously one-HSP hit can't have the evidence we look for
        return False
    valid_hsps = filter(lambda x: abs(x.hit_pos[1]-x.hit_pos[0]) > len_cutoff,
                        hit.hsps)
    for h1, h2 in combinations(valid_hsps, 2):
        if overlap(h1.hit_pos, h2.hit_pos) and \
                not overlap(h1.query_pos, h2.query_pos):
            l = overlap_len(h1.hit_pos, h2.hit_pos)
            if l >= abs(h1.hit_pos[1]-h1.hit_pos[0]) * overlap_cutoff and \
                    l >= abs(h2.hit_pos[1]-h2.hit_pos[0]) * overlap_cutoff:
                return True
    return False
