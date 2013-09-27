
__author__ = 'Soh Ishiguro <t10078si@sfc.keio.ac.jp>'
__version__ = '0.0.1'

import string
import os.path
from Bio import SeqIO


class SequenceIOException(Exception):
    def __init__(self, file):
        self.file = file
        
    def __str__(self):
        return ""


class SequenceIO(object):
    def __init__(self, file='', circular=True):

        if os.path.splitext(file)[1] != '.fasta':
            raise RuntimeError("Error: fasta format is only accepted")
        if not os.path.exists(file):
            raise RuntimeError("Error: given file is not found")
            
        self.seq = file
        self.circular = circular # default: circular genome
    
    def read_sequence(self):
        record = SeqIO.read(self.seq, "fasta")
        self.rec = []
        for i in record.seq:
            self.rec.append([i])
        return self.rec
    
    def insert_sequence(self, start, original_seq, insert_seq):
        pre = original_seq[0:start]
        post = original_seq[start:end]
        self.inesrt = pre + insert_seq + post
        return self.insert
        
    def sequence_file_name(self):
        return self.seq
    
    def fasta_header(self):
        record = SeqIO.read(self.seq, "fasta")
        return record.id

    def get_complement_seq(self):
        self.rev_seq = []
        record = SeqIO.read(self.seq, 'fasta')
        
        for nuc in record.seq:
            if nuc == 'A' or nuc == 'a':
                self.rev_seq.append(['T'])
            elif nuc == 'T' or nuc == 't':
                self.rev_seq.append(['A'])
            elif nuc == 'G' or nuc == 'g':
                self.rev_seq.append(['C'])
            elif nuc == 'C' or nuc == 'c':
                self.rev_seq.append(['G'])
            elif nuc == 'n' or nuc == 'N':
                self.rev_seq.append(['n'])
                
        return self.rev_seq
        
    def get_reverse_complement_seq(self):
        self.rev_comp_seq = []
        record = SeqIO.read(self.seq, 'fasta')
        return self.rev_comp_seq

    def print_sequence(self, all=False):
        limit = 100
        concentrated_seq = ''
        
        if all:
            for s in self.rec:
                concentrated_seq += str(s[0])
            return concentrated_seq
   
        elif not all:
            for index, seq in enumerate(self.rec):
                concentrated_seq += str(seq[0])
                if index > limit:
                    break
            return concentrated_seq
        


