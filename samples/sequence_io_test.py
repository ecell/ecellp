#!/usr/bin/python 

from sequence_io import SequenceIO
from genomic_coordinate import GenomicCoordinate
from genomic_attribute import GenomicAttribute

def test_sequence_io_class():
    seq_file = '../../data/seq.fasta'
    seq_io = SequenceIO(file=seq_file, circular=True)
    print "fasta_header() is called, %s" %(seq_io.fasta_header())

    genomic_array = seq_io.read_sequence()
    print "read_sequence() is called"

    print "print_sequence() is called, %s" %(seq_io.print_sequence(all=False))
    print "sequence_file_name() is called, %s" %(seq_io.sequence_file_name())
    print "get_complement_seq() is called, %s" %(seq_io.get_complement_seq())
    print "get_reverse_complement_seq() is called, %s" %(seq_io.get_reverse_complement_seq())

def test_genomic_coordinate_class():
    pass

def test_genomic_attribute_class():
    pass


if __name__ == '__main__':
    test_sequence_io_class()





