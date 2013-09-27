
__author__ = 'Soh Ishiguro <t10078si@sfc.keio.ac.jp>'
__version__ = '0.0.1'


class GenomicAttribute(object):
    def __init__(self, start='', end=''):
        self.start = start
        self.end = end

    def get_gene_seq(self):
        pass

    def get_gene_name(self):
        pass

    def get_gene_annotation_from_position(self):
        pass

    def find_all_gene_names(self):
        pass

    def find_all_cds(self):
        pass

    def is_strand(self):
        pass

    def is_circular(self):
        pass

    def size_of_genome(self):
        pass

    def is_traslated(self):
        pass


