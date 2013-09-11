__program__ = 'core'
__version__ = '0.1'
__author__ = 'Kazunari Kaizu <kaizu@riken.jp>'
__copyright__ = ''
__license__ = ''

import copy
import numpy


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

def get_region(seq, *region):
    (start, end, stride, is_cyclic) = region
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

def set_region(seq, value, *region):
    map_region(lambda x: value, seq, *region)

def map_region(func, seq, *region):
    (start, end, stride, is_cyclic) = region
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

class Alignment(object):

    def __init__(self, seq, is_cyclic = True, **annotations):
        self.__seq = seq
        self.__is_cyclic = is_cyclic
        self.set_annotations(**annotations)

        self.__attributes = [[] for _ in range(len(self.__seq))]

    def is_cyclic(self):
        return self.__is_cyclic

    def annotations(self):
        return self.__annotations

    def set_annotations(self, **kwargs):
        self.__annotations = kwargs

    def size(self):
        return len(self.__seq)

    def num_domains(self):
        return len(self.list_domains())

    def list_domains(self, dom = None):
        if dom is None:
            return list(set(sum(self.__attributes, [])))
        else:
            return [elem for elem in set(sum(self.__attributes, []))
                    if dom.match(elem, self)]
            # return [elem for elem in set(sum(self.__attributes, []))
            #         if elem.name() == dom.name()]

    def place_domain(self, dom, *coord):
        value = copy.deepcopy(dom)
        (start, end, stride) = coord
        value.set_attributes(start=start, end=end, stride=stride)
        func = lambda attrs: attrs + [value]
        map_region(func, self.__attributes, start, end, stride, self.is_cyclic())

    def query_domains_by_region(self, *coord):
        (start, end, stride) = coord
        domains = get_region(self.__attributes, start, end, stride, self.is_cyclic())
        return list(set(sum(domains, [])))

    def sequence(self, *region):
        if len(region) == 0:
            return self.__seq
        elif len(region) == 3:
            (start, end, stride) = region
            if stride > 0:
                return get_region(
                    self.__seq, start, end, stride, self.is_cyclic())
            else:
                return reverse_sequence(
                    get_region(
                        self.__seq, start, end, stride, self.is_cyclic()))
        else:
            raise ValueError, "Invalid arguments [%s]" % (str(region))

    def find(self, seq, stride):
        if stride > 0:
            idx = self.sequence().find(seq)
            if idx < 0:
                return None
            else:
                return (idx + 1, idx + len(seq), stride)
        elif stride < 0:
            idx = reverse_sequence(self.sequence()).find(seq)
            if idx < 0:
                return None
            else:
                return (self.size() - idx, self.size() - (idx + len(seq)) + 1, stride)
        else:
            raise ValueError, "Invalid argument given [%s]" % (stride)

class Track(object):

    def __init__(self, size, is_cyclic = True):
        self.__size = size
        self.__is_cyclic = is_cyclic

        self.__states = numpy.zeros(size, numpy.uint64)
        self.__idgen = IDGenerator()
        self.__domains = {}
        self.__dserial_did_map = {}

    def size(self):
        return self.__size

    def is_cyclic(self):
        return self.__is_cyclic

    def states(self):
        return self.__states

    def list_domains(self, dom = None):
        retval = []
        if dom is None:
            for i in range(self.size()):
                did = self.__states[i]
                if did == 0:
                    continue
                elem = copy.deepcopy(self.__domains[did])
                elem.set_attributes(start = i + 1, end = i + 1)
                retval.append(elem)
        else:
            did = self.__dserial_did_map.get(dom.serial())
            if did is None:
                return retval
            for i in range(self.size()):
                if self.__states[i] != did:
                    continue
                elem = copy.deepcopy(self.__domains[did])
                elem.set_attributes(start = i + 1, end = i + 1)
                retval.append(elem)
        return retval

    def num_molecules(self, dom):
        did = self.__dserial_did_map.get(dom.serial())
        if did is None:
            return 0
        else:
            return sum(self.__states == did)

    def place_molecule(self, dom, coord):
        self.place_domain(dom, coord)

    def remove_molecule(self, coord):
        return self.remove_domain(coord)[0]

    def place_domain(self, dom, *coord):
        if len(coord) == 3:
            region = (coord[0], coord[1], coord[2], self.is_cyclic())
        elif len(coord) == 1:
            region = (coord[0], coord[0], +1, self.is_cyclic())
        else:
            raise RuntimeError, "Invalid arguments given [%s]" % (str(coord))

        retval = get_region(self.__states, *region)
        if any(retval != 0):
            raise RuntimeError, "Already occupied."

        did = self.__dserial_did_map.get(dom.serial())
        if did is None:
            did = self.__idgen()
            self.__domains[did] = dom
            self.__dserial_did_map[dom.serial()] = did

        set_region(self.__states, did, *region)

    def remove_domain(self, *coord):
        if len(coord) == 3:
            region = (coord[0], coord[1], coord[2], self.is_cyclic())
        elif len(coord) == 1:
            region = (coord[0], coord[0], +1, self.is_cyclic())
        else:
            raise RuntimeError, "Invalid arguments given [%s]" % (str(coord))

        substates = get_region(self.__states, *region).copy()
        if any(substates == 0):
            raise RuntimeError, "No domain assigned."

        set_region(self.__states, 0, *region)

        retval = []
        for did in set(substates):
            if all(self.__states != did):
                dom = self.__domains.pop(did)
                del self.__dserial_did_map[dom.serial()]
                retval.append(dom)
            else:
                retval.append(self.__domains[did])

        return retval

        # pos = resolve_coordinate(coord, self.size(), self.is_cyclic())
        # did = self.__states[pos]
        # if did == 0:
        #     raise RuntimeError, "No domain assigned."
        # self.__states[pos] = 0

        # if any(self.__states == did):
        #     return self.__domains[did]
        # else:
        #     dom = self.__domains.pop(did)
        #     del self.__dserial_did_map[dom.serial()]
        #     return dom

    # def remove_domain(self, dom):
    #     did = self.__dserial_did_map.get(dom.serial())
    #     if did is None:
    #         return

    #     for i in range(len(self.__states)):
    #         if self.__states[i] == did:
    #             self.__states[i] = 0

    #     dom = self.__domains.pop(did)
    #     del self.__dserial_did_map[dom.serial()]

    # def get_domain(self, coord):
    #     pos = resolve_coordinate(coord, self.size(), self.is_cyclic())
    #     did = self.__states[pos]
    #     if did == 0:
    #         return None
    #     else:
    #         return self.__domains[did]

    def query_domains_by_region(self, *coord):
        (start, end, stride) = coord
        dids = get_region(self.__states, start, end, stride, self.is_cyclic())
        return list(set([self.__domains[did] for did in dids if did != 0]))

class AlignmentSpace(object):

    def __init__(self, algn):
        self.__alignment = algn
        self.__tracks = []
        self.__tname_idx_map = {}

    def size(self):
        return self.__alignment.size()

    def is_cyclic(self):
        return self.__alignment.is_cyclic()

    def alignment(self):
        return self.__alignment

    def create_track(self, name):
        if name in self.__tname_idx_map.keys():
            raise RuntimeError, "Track [%s] already exists" % name
        trk = Track(self.size(), self.is_cyclic())
        self.__tname_idx_map[name] = len(self.__tracks)
        self.__tracks.append(trk)
        return trk

    def list_domains(self, *args, **kwargs):
        retval = self.__alignment.list_domains(*args, **kwargs)
        for trk in self.__tracks:
            retval.extend(trk.list_domains(*args, **kwargs))
        return retval

    def query_domains_by_region(self, *args, **kwargs):
        retval = self.__alignment.query_domains_by_region(*args, **kwargs)
        for trk in self.__tracks:
            retval.extend(trk.query_domains_by_region(*args, **kwargs))
        return retval

class CompartmentSpace(object):

    def __init__(self, volume):
        self.__volume = volume

        self.__domains = []
        self.__num_molecules = []
        self.__dserial_idx_map = {}

    def volume(self):
        return self.__volume

    def list_domains(self):
        return self.__domains

    def num_molecules(self, dom):
        idx = self.__dserial_idx_map.get(dom.serial())
        if idx is None:
            return 0
        else:
            return self.__num_molecules[idx]

    def add_molecules(self, dom, num):
        idx = self.__dserial_idx_map.get(dom.serial())
        if idx is None:
            idx = len(self.__domains)
            self.__dserial_idx_map[dom.serial()] = idx
            self.__domains.append(dom)
            self.__num_molecules.append(0)

        self.__num_molecules[idx] += num

    def remove_molecules(self, dom, num):
        idx = self.__dserial_idx_map.get(dom.serial())
        if idx is None:
            raise RuntimeError, "Domain [%s] not found" % (str(dom))
        elif self.__num_molecules[idx] < num:
            raise RuntimeError, ""

        self.__num_molecules[idx] -= num
        if self.__num_molecules[idx] == 0:
            if idx != len(self.__domains) - 1:
                tmp = self.__domains[-1]
                self.__dserial_idx_map[tmp.serial()] = idx
                self.__domains[idx] = tmp
                self.__num_molecules[idx] = self.__num_molecules[-1]
            self.__domains.pop()
            self.__num_molecules.pop()
            del self.__dserial_idx_map[dom.serial()]

class IDGenerator(object):

    def __init__(self, stride = 0):
        self.__stride = stride

    def __call__(self):
        self.__stride += 1
        return self.__stride

class Domain(object):

    def __init__(self, name, **attrs):
        self.__name = name
        self.__attributes = {}

        self.set_attributes(**attrs)

    def set_attributes(self, **attrs):
        if "name" in attrs.keys():
            raise RuntimeError, 'Do not overload "name"'
        self.__attributes.update(attrs)

    def get_attribute(self, key):
        return self.__attributes.get(key)

    def attributes(self):
        return self.__attributes

    def serial(self):
        return self.name()

    def name(self):
        return self.__name

    def match(self, dom, parent):
        return dom.name() == self.name()

    def as_Domain(self):
        return Domain(self.name(), **self.attributes())

    def __repr__(self):
        mod = self.__class__.__module__
        cls = self.__class__.__name__
        # mem = '0x' + hex(id(self))[2:].zfill(8).upper()
        return '<{0}.{1} instance with serial "{2}">'.format(mod, cls, self.serial())

class Condition(object):

    def __init__(self, left, right, includes = [], excludes = []):
        if len(includes) == 0 and len(excludes) == 0:
            raise RuntimeError, ""

        self.__bounds = (left, right)
        self.__includes = includes
        self.__excludes = excludes

    def match(self, dom, parent):
        start = dom.get_attribute("start")
        stride = dom.get_attribute("stride")
        region = (
            start + self.__bounds[0] * stride,
            start + self.__bounds[1] * stride,
            stride)
        neighbors = parent.query_domains_by_region(*region)

        # do more precise filtering before checking

        for elem in self.__includes:
            if not any([elem.match(neighbor, parent) for neighbor in neighbors]):
                return False
        for elem in self.__excludes:
            if any([elem.match(neighbor, parent) for neighbor in neighbors]):
                return False
        return True

class FilteringDomain(Domain):

    def __init__(self, name = "", **attrs):
        Domain.__init__(self, name, **attrs)

    def match(self, dom, parent):
        if self.name() is not "" and self.name() != dom.name():
            return False

        for key, value in self.attributes().items():
            if dom.get_attribute(key) != value:
                return False
        else:
            return True

class ConditionalDomain(Domain):

    def __init__(self, name, **attrs):
        Domain.__init__(self, name, **attrs)

        self.__conditions = []

    def add_condition(self, condition):
        self.__conditions.append(condition)

    def find(self, parent):
        retval = []
        for dom in parent.list_domains(self.as_Domain()):
            if self.match(dom, parent):
                retval.append(dom)
        return retval

    def match(self, dom, parent):
        if not self.as_Domain().match(dom, parent):
            return False

        for condition in self.__conditions:
            if not condition.match(dom, parent):
                return False
        else:
            return True

class SecondOrderInteraction(object):

    def __init__(self, r1, r2, p1, k):
        self.__reactants = (r1, r2)
        self.__products = (p1, )
        self.__k = k

    def reactants(self):
        return self.__reactants

    def products(self):
        return self.__products

    def k(self):
        return self.__k


if __name__ == "__main__":
    pass
