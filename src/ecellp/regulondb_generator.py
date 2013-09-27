import re
import os.path

import utils
import species


class RegulonDBPromoterDecGenerator(object):

    def __init__(self, filename):
        self.__url = ("http://regulondb.ccg.unam.mx/"
                      "menu/download/datasets/files/PromoterSet.txt")
        self.__filename = filename

    def setup(self):
        if os.path.isfile(self.__filename):
            return True
        else:
            return utils.fetch_url(self.__url, self.__filename)

    def reformat(self, data):
        if data[2] in "forward":
            return (data[0], "+1", data[3], data[3], "promoter", data[5])
        else:
            return (data[0], "-1", data[3], data[3], "promoter", data[5])

    def generate(self, session):
        if not self.setup():
            raise RuntimeError, "File [%s] not found." % (self.__filename)

        with open(self.__filename, "r") as fin:
            while True:
                line = fin.readline()
                if line is None or line == "":
                    break
                elif line.isspace() or re.match('^#', line):
                    continue

                line = line.rstrip()

                data = line.split("\t")
                if not len(data) == 7:
                    continue

                # fout.write("%s\n" % ("\t".join(self.reformat(data))))
                obj = species.PromoterDec(*self.reformat(data))
                session.add(obj)
        session.commit()

class RegulonDBTerminatorDecGenerator(object):

    def __init__(self, filename):
        self.__url = ("http://regulondb.ccg.unam.mx/"
            "menu/download/datasets/files/TerminatorSet.txt")
        self.__filename = filename

    def setup(self):
        if os.path.isfile(self.__filename):
            return True
        else:
            return utils.fetch_url(self.__url, self.__filename)

    def reformat(self, data):
        if data[3] in "forward":
            return (data[0], "+1", data[1], data[2], "terminator", data[4])
        else:
            return (data[0], "-1", data[1], data[2], "terminator", data[4])

    def generate(self, session):
        if not self.setup():
            raise RuntimeError, "File [%s] not found." % (self.__filename)

        with open(self.__input_filename, "r") as fin:
            while True:
                line = fin.readline()
                if line is None or line == "":
                    break
                elif line.isspace() or re.match('^#', line):
                    continue

                line = line.rstrip()

                if not re.match('^\S',line):
                    continue

                data = line.split("\t")
                # fout.write("%s\n" % ("\t".join(self.reformat(data))))
                obj = species.TerminatorDec(*self.reformat(data))
                session.add(obj)
        session.commit()

class RegulonDBOperonDecGenerator(object):

    def __init__(self, filename):
        filename = 'data/OperonSet.txt'
        self.__filename = filename

    def reformat(self, data):
        name = data[0]
        start = int(data[1])
        end  = int(data[2])
        strand = data[3]
        num_operon = int(data[4])
        gene_name = data[5] #XXX
        sequence = "" #XXX

        if 'forward' in strand:
            return (name, +1, start, end, "operon", sequence)
        elif 'reverse' in strand:
            return (name, -1, start, end, "operon", sequence)

    def generate(self, session):
        """This script generates the operon DB from OperonSet.txt,
        input data was collected from web resource of the RegulonDB.
        URL    : http://regulondb.ccg.unam.mx/menu/download/datasets/files/OperonSet.txt
        OUTPUT : promoter_name    start    end    strand    gene_name(s)
        """
        header = "promoter_id\tstrand\tstart\tend\tfeature\tsequence"

        with open(self.__filename, 'r') as fin:
            for line in fin:
                if line.startswith("#") or line.isspace():
                    continue

                data = line.strip().split("\t")
                obj = species.OperonDec(*self.reformat(data))
                session.add(obj)
        session.commit()

class RegulonDBPromoterInteractionDecGenerator(object):
    """This class should be integrated with RegulonDBOperonDecGenerator
    in future."""

    def __init__(self, filename):
        filename = 'data/OperonSet.txt'
        self.__filename = filename

    def reformat(self, data):
        name = data[0]
        start = int(data[1])
        end  = int(data[2])
        strand = data[3]
        num_operon = int(data[4])
        gene_name = data[5] #XXX
        sequence = "" #XXX

        if 'forward' in strand:
            return (name, +1, start, end, "operon", sequence)
        elif 'reverse' in strand:
            return (name, -1, start, end, "operon", sequence)

    def generate(self, session):
        """This script is used for generating promoter information per gene
        DATA URL: http://regulondb.ccg.unam.mx/menu/download/datasets/files/OperonSet.txt
        OUTPUT  : gene_name    strand    operon_id
        """
        gene_map = {}
        with open(self.__filename, 'r') as fin:
            for line in fin:
                if line.startswith("#") or line.isspace():
                    continue

                data = line.strip().split("\t")
                args = self.reformat(data)

                operon_name = 'OPERON_' + args[0]
                strand = args[1]
                gene_names  = data[5]

                genes = gene_names.split(",")
                for gene in genes:
                    gene_map.update(
                        {gene : {'operon_id' : operon_name, 'strand' : strand}})

        raise NotImplementedError, "not fully implemented yet."
