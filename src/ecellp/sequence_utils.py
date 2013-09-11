__program__ = 'sequence_utils'

__version__ = '0.1'
__author__ = 'Kazunari Kaizu <kaizu@riken.jp>'
__copyright__ = ''
__license__ = ''


def read_sequence(filename):
    seq = ""
    with open(filename, "r") as fin:
        while True:
            line = fin.readline()
            if line == "":
                break
            line = line.strip()
            seq += line
    return seq

def format_sequence(seq):
    formatted = ""
    start, stride1, stride2 = 0, 10, 6
    while start < len(seq):
        s = ' '.join([
            seq[start + i * stride1: start + (i + 1) * stride1]
            for i in range(stride2)])
        s += " " * ((stride1 + 1) * stride2 - 1 - len(s))
        formatted += "%05d %s %05d\n" % (start + 1, s, start + stride1 * stride2)
        start += stride1 * stride2
    return formatted[: -1]

def resolve_coordinate(coord, len_seq, is_cyclic = False):
    if not is_cyclic:
        if not 0 < coord <= len_seq:
            raise ValueError, "Invalid coordinate was given [%s]." % (coord)
    elif coord == 0:
        raise ValueError, "A coordinate must be non-zero [%s]." % (coord)
    return (coord - 1) % len_seq

def get_region(seq, region, is_cyclic = False):
    (start, end, stride) = region
    start, end = (
        resolve_coordinate(start, len(seq), is_cyclic),
        resolve_coordinate(end, len(seq), is_cyclic))

    if stride > 0:
        if end >= start:
            return seq[start: end + 1: stride]
        elif is_cyclic:
            return seq[start: : stride] + seq[: end + 1: stride]
        else:
            raise RuntimeError, "Invalid coordinate was given [%s]" % (
                region[: -1], )
    elif stride < 0:
        if start >= end:
            return seq[end: start + 1: -stride]
        elif is_cyclic:
            return seq[end: : -stride] + seq[: start + 1: -stride]
        else:
            raise RuntimeError, "Invalid coordinate was given [%s]" % (
                region[: -1], )
    else:
        raise RuntimeError, "A stride must be non-zero."

def set_region(seq, value, region, is_cyclic = False):
    map_region(lambda x: value, seq, region, is_cyclic)

def map_region(func, seq, region, is_cyclic = False):
    (start, end, stride) = region
    start, end = (
        resolve_coordinate(start, len(seq), is_cyclic),
        resolve_coordinate(end, len(seq), is_cyclic))

    if stride > 0:
        if end >= start:
            for i in range(start, end + 1, stride):
                seq[i] = func(seq[i])
        elif is_cyclic:
            for i in range(start, len(seq), stride) + range(0, end + 1, stride):
                seq[i] = func(seq[i])
        else:
            raise RuntimeError, "Invalid coordinate was given [%s]" % (
                region[: -1], )
    elif stride < 0:
        if start >= end:
            for i in range(end, start + 1, -stride):
                seq[i] = func(seq[i])
        elif is_cyclic:
            for i in range(end, len(seq), -stride) + range(0, start + 1, -stride):
                seq[i] = func(seq[i])
        else:
            raise RuntimeError, "Invalid coordinate was given [%s]" % (
                region[: -1], )
    else:
        raise RuntimeError, "A stride must be non-zero."

def reverse_sequence(seq):
    return ''.join([
        'A' if c is 'T' else 'T' if c is 'A' else 'G' if c is 'C' else 'C'
        for c in seq[: : -1]])
