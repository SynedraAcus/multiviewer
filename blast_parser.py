"""
A lightweight BLAST tabular output parser.
It slightly outperforms Biopython's BLAST parser, but is mostly here as an
exercise in a more functional programming style and proper unittests.
"""

from collections import namedtuple
from io import TextIOBase


BlastHSP = namedtuple('BlastHSP', ['query_id', 'hit_id', 'query_pos', 'hit_pos',
                                   'evalue'])
BlastHit = namedtuple('BlastHit', ['query_id', 'hit_id', 'hsps'])


def create_handle(filename):
    """
    Take a readable filehandle or an str, return readable filehandle.
    If it's a filehandle already, this is an identity function, otherwise the
    string is treated as a filename to be open. Any other argument type raises
    TypeError.
    :param filename:
    :return:
    """
    if isinstance(filename, str):
        return open(filename, mode='r')
    elif isinstance(filename, TextIOBase) and filename.readable():
        return filename
    else:
        raise TypeError('Only string or readable text-mode filehandle accepted by parse_blast_file')
    
    
def parse_blast_line(line):
    """
    Take a single line from BLAST file and return a BlastHSP.
    It assumes the default tabular output from BLAST 2.2.28+ (generated with
    `-outfmt 6` without any additional format specifications). May or may not
    work with other BLAST versions.
    :param line: str
    :return: BlastHSP
    """
    # Last field is ignored anyway
    # line = line.rstrip('\n')
    arr = line.split('\t')
    if not len(arr) == 12:
        raise ValueError('Incorrect BLAST line {}'.format(line))
    #  There is no explicit validity check, but hopefully IndexError or some
    #  exception of the type converters will be raised on the invalid string
    return BlastHSP(query_id=arr[0],
                    hit_id=arr[1],
                    query_pos=(int(arr[6]), int(arr[7])),
                    hit_pos=(int(arr[8]), int(arr[9])),
                    evalue=float(arr[10]))


def assemble_hits(iterable):
    """
    Groups HSPs into BlastHit instances
    Takes an iterable of BlastHSPs, yields BlastHits. Assumes the iterable to be
    sorted by query ID and, for each query, to be sorted by hit ID (*vice versa
    should also work).
    :param iterable:
    :return:
    """
    y = []
    current_hit = ''
    current_query = ''
    for blast_hsp in iterable:
        #  The following two checks will only run when the generator is created
        if not current_hit:
            current_hit = blast_hsp.hit_id
        if not current_query:
            current_query = blast_hsp.query_id
        if current_hit == blast_hsp.hit_id and current_query == blast_hsp.query_id:
            y.append(blast_hsp)
        else:
            yield BlastHit(query_id=current_query, hit_id=current_hit, hsps=y)
            current_hit = blast_hsp.hit_id
            current_query = blast_hsp.query_id
            y = [blast_hsp]
    yield BlastHit(query_id=current_query, hit_id=current_hit, hsps=y)


def parse_blast_file_to_hsps(filename, ignore_trivial=True):
    """
    Take a BLAST output file and generate BlastHSP instances from its lines.
    It assumes the default tabular output from BLAST 2.2.28+ (generated with
    `outfmt 6` or `outfmt 7` without any additional format specifications).
    Any line starting with `#` is considered a comment and omitted.
    `file` can be either a filehandle in text mode or a filename. In latter case
    it would be opened in `'r'` mode.
    By default ignores trivial hits (ie ones where hit and query names are
    identical).
    :param filename: str or filehandle
    :param ignore_trivial: bool
    :return:
    """
    for line in create_handle(filename):
        if not line[0] == '#':
            hsp = parse_blast_line(line)
            if (not ignore_trivial) or (not hsp.query_id == hsp.hit_id):
                yield hsp
    # Not closing a filehandle because it may be used by the calling code.


def parse_blast_file_to_hits(filename, ignore_trivial=True):
    """
    Take a BLAST output file and generate BlastHit instances from its lines.
    It assumes the default tabular output from BLAST 2.2.28+ (generated with
    `outfmt 6` or `outfmt 7` without any additional format specifications).
    Any line starting with `#` is considered a comment and omitted.
    `file` can be either a filehandle in text mode or a filename. In latter case
    it would be opened in `'r'` mode.
    :param filename: str or filehandle
    :param ignore_trivial: bool
    :return: generator
    """
    return assemble_hits(parse_blast_file_to_hsps(filename, ignore_trivial))


def iterate_by_query(hits_iterator):
    """
    Yield lists of hits with the same query.
    Given an iterator of BlastHit instances, this function produces lists of the
    BlastHits with the same query. It assumes hits to be sorted by query, ie
    that all the hits of a given query are following one after another in the
    list (which is a normal behaviour for BLAST output)
    :param hits_iterator:
    :return:
    """
    current_query = ''
    #  Probably could be an iterator as well, but I can't design the algorithm.
    hits = []
    for blasthit in hits_iterator:
        if not current_query:
            current_query = blasthit.query_id
        if blasthit.query_id == current_query:
            hits.append(blasthit)
        else:
            yield hits
            hits = [blasthit]
            current_query = blasthit.query_id
    yield hits
