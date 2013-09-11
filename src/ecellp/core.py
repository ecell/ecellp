__program__ = 'core'
__version__ = '0.1'
__author__ = 'Kazunari Kaizu <kaizu@riken.jp>'
__copyright__ = ''
__license__ = ''

import copy
import numpy

from sequence_utils import *


class IDGenerator(object):

    def __init__(self, stride = 0):
        self.__stride = stride

    def __call__(self):
        self.__stride += 1
        return self.__stride

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
        map_region(func, self.__attributes, coord, self.is_cyclic())

    def query_domains_by_region(self, *region):
        domains = get_region(self.__attributes, region, self.is_cyclic())
        return list(set(sum(domains, [])))

    def sequence(self, *region):
        if len(region) == 0:
            return self.__seq
        elif len(region) == 3:
            (_, _, stride) = region
            if stride > 0:
                return get_region(
                    self.__seq, region, self.is_cyclic())
            else:
                return reverse_sequence(
                    get_region(
                        self.__seq, region, self.is_cyclic()))
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

    # def num_molecules(self, dom):
    #     did = self.__dserial_did_map.get(dom.serial())
    #     if did is None:
    #         return 0
    #     else:
    #         return sum(self.__states == did)

    # def place_molecule(self, dom, coord):
    #     self.place_domain(dom, coord)

    # def remove_molecule(self, coord):
    #     return self.remove_domain(coord)[0]

    def place_domain(self, dom, *coord):
        if len(coord) == 3:
            region = coord
        elif len(coord) == 1:
            region = (coord[0], coord[0], +1)
        else:
            raise RuntimeError, "Invalid arguments given [%s]" % (str(coord))

        retval = get_region(self.__states, region, self.is_cyclic())
        if any(retval != 0):
            raise RuntimeError, "Already occupied."

        did = self.__dserial_did_map.get(dom.serial())
        if did is None:
            did = self.__idgen()
            self.__domains[did] = dom
            self.__dserial_did_map[dom.serial()] = did

        set_region(self.__states, did, region, self.is_cyclic())

    def remove_domain(self, *coord):
        if len(coord) == 3:
            region = coord
        elif len(coord) == 1:
            region = (coord[0], coord[0], +1)
        else:
            raise RuntimeError, "Invalid arguments given [%s]" % (str(coord))

        substates = get_region(self.__states, region, self.is_cyclic()).copy()
        if any(substates == 0):
            raise RuntimeError, "No domain assigned."

        set_region(self.__states, 0, region, self.is_cyclic())

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
        substates = get_region(self.__states, coord, self.is_cyclic())
        return list(set([self.__domains[did] for did in substates if did != 0]))

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
        retval = 0
        for elem, num in zip(self.__domains, self.__num_molecules):
            if dom.match(elem):
                retval += num
        return retval

        # idx = self.__dserial_idx_map.get(dom.serial())
        # if idx is None:
        #     return 0
        # else:
        #     return self.__num_molecules[idx]

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
