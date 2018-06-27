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
        # TODO: string export
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

