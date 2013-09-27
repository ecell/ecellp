
__author__ = 'Soh Ishiguro <t10078si@sfc.keio.ac.jp>'
__version__ = '0.0.1'

import string
import os.path


class GenomicCoordinateException(Exception):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __str__(self):
        return "Genomic position is only 1-origin, given value of %d, %d were not accepted" % (self.start, self.end)


class GenomicCoordinate(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def is_circular(self):
        return self.start
    
    def seq_name(self):
        return self.fasta

    def get_sequence_region_from_circular(self):
        if self.start > self.end:
            return self.sequence[self.start-1:self.end]

        elif self.end < self.start:
            return self.sequence[self.start:self.length % self.end]
            
        elif self.start > 0 and self.end > 0:
            return self.sequence[self.start-1:self.end]

        elif self.start > 0 and self.end < 0:
            return self.sequence[self.start-1:self.length-abs(self.end)]

        elif self.start < 0 and self.end > 0:
            return self.sequence[self.start-1:self.length-abs(self.end)]
            
        elif self.start < 0 and self.end < 0:
            return self.sequence[self.length-abs(end):self.length-abs(self.start)]
        
    def get_sequence_region_from_linear(self):
        if self.start > 0 and self.end > 0:
            return self.sequence[self.start-1:self.end]

        elif self.start > 0 and self.end > self.length:
            return self.sequence[self.start-1:self.length]
            
        elif self.start > self.length:
            return ''

    def map_attribute_to_genome(self):
        pass



