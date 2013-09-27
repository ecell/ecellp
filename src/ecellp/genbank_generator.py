import re
import os

from Bio import SeqIO

import species


class GenbankDecGenerator(object):

    def __init__(self, filename):
        self.__filename = filename

    def generate(self, session):
        with open(self.__filename, "r") as fin:
            for record in SeqIO.parse(fin, 'genbank'):
                for feature in record.features:
                    if feature.type == 'CDS':
                        if (not 'pseudo' in feature.qualifiers
                            and 'gene' in feature.qualifiers):
                            name = feature.qualifiers['gene'][0]
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            strand = int(feature.location.strand)
                            seq = str(record.seq[start: end])
                            feature = feature.type

                            data = (name, strand, start, end, feature, seq)
                            session.add(species.CDSDec(*data))
                    elif feature.type == 'rRNA':
                        name = feature.qualifiers['gene'][0]
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = int(feature.location.strand)
                        seq = str(record.seq[start: end])
                        feature = feature.type

                        data = (name, strand, start, end, feature, seq)
                        session.add(species.rRNADec(*data))
                    elif feature.type == 'tRNA':
                        name = feature.qualifiers['gene'][0]
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = int(feature.location.strand)
                        seq = str(record.seq[start: end])
                        feature = feature.type

                        data = (name, strand, start, end, feature, seq)
                        session.add(species.tRNADec(*data))
        session.commit()
