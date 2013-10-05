__program__ = 'session'
__version__ = '0.0.1'
__author__ = 'Soh Ishiguro <t10078si@sfc.keio.ac.jp>'
__copyright__ = ''
__license__ = ''

import sys
import os
import ConfigParser

from sqlalchemy import *
from sqlalchemy.orm import sessionmaker

import species
import genbank_generator
import regulondb_generator


def generate_decs(session, conf):
    gen = genbank_generator.GenbankDecGenerator(conf.GENBANK_FILE)
    gen.generate(session)

    gen = regulondb_generator.RegulonDBPromoterDecGenerator(
        conf.PROMOTER_FILE)
    gen.generate(session)

    gen = regulondb_generator.RegulonDBPromoterDecGenerator(
        conf.TERMINATOR_FILE)
    gen.generate(session)

path = os.path.abspath(os.path.dirname(__file__))
sys.path.append(path)

class DBConfig(object):

    def __init__(self, filename=None, paths=None):
        '''
        Arguments:
            filename: OPTIONAL, indicates conf.ini path.
            paths:    OPTIONAL, a dictionary which contains paths to data.
        '''

        APP_ROOT = self.find_app_root()

        # Input files
        self.GENBANK_FILE    = "/data/NC_000913.gbk"
        self.PROMOTER_FILE   = "/data/PromoterSet.txt"
        self.TERMINATOR_FILE = "/data/TerminatorSet.txt"
        self.GENOME_SEQUENCE = "/data/test.fa"

        # Output files
        # self.CDS_OUT         = "/data/CDS_annotation.tbl"
        # self.rRNA_OUT        = "/data/rRNA_annotation.tbl"
        # self.tRNA_OUT        = "/data/tRNA_annotation.tbl"
        # self.PROMOTER_OUT    = "/data/promoter_annotation.tbl"
        # self.TERMINATOR_OUT  = "/data/terminator_annotation.tbl"

        # DB
        self.DB_PATH = "/data/ecoli.sqlite3"

        if filename:
            self.__filename = filename
            if not os.path.isfile(self.__filename):
                raise RuntimeError, "Configuration file [%s] is not found" % (
                    os.path.abspath(self.__filename))
            conf = ConfigParser.RawConfigParser()
            conf.read(self.__filename)

            if conf.has_option('input_data', 'genbank'):
                self.GENBANK_FILE    = conf.get('input_data', 'genbank')
            if conf.has_option('input_data', 'promoter'):
                self.PROMOTER_FILE   = conf.get('input_data', 'promoter')
            if conf.has_option('input_data', 'terminator'):
                self.TERMINATOR_FILE = conf.get('input_data', 'terminator')
            if conf.has_option('input_data', 'sequence'):
                self.GENOME_SEQUENCE = conf.get('input_data', 'sequence')

            # if conf.has_option('output_data', 'cds'):
            #     self.CDS_OUT         = conf.get('output_data', 'cds')
            # if conf.has_option('output_data', 'rrna'):
            #     self.rRNA_OUT        = conf.get('output_data', 'rrna')
            # if conf.has_option('output_data', 'trna'):
            #     self.tRNA_OUT        = conf.get('output_data', 'trna')
            # if conf.has_option('output_data', 'promoter'):
            #     self.PROMOTER_OUT    = conf.get('output_data', 'promoter')
            # if conf.has_option('output_data', 'terminator'):
            #     self.TERMINATOR_OUT  = conf.get('output_data', 'terminator')

            if conf.has_option('db', 'db_path'):
                self.DB_PATH         = conf.get('db', 'db_path')

        if paths:
            if paths.has_key('genbank'):
                self.GENBANK_FILE    = paths['genbank']
            if paths.has_key('promoter_file'):
                self.PROMOTER_FILE   = paths['promoter_file']
            if paths.has_key('terminator_file'):
                self.TERMINATOR_FILE = paths['terminator_file']
            if paths.has_key('sequence'):
                self.GENOME_SEQUENCE = paths['sequence']

            # if paths.has_key('cds_out'):
            #     self.CDS_OUT         = paths['cds_out']
            # if paths.has_key('rrna_out'):
            #     self.rRNA_OUT        = paths['rrna_out']
            # if paths.has_key('trna_out'):
            #     self.tRNA_OUT        = paths['trna_out']
            # if paths.has_key('promoter_out'):
            #     self.PROMOTER_OUT    = paths['promoter_out']
            # if paths.has_key('terminator_out'):
            #     self.TERMINATOR_OUT  = paths['terminator_out']

            if paths.has_key('db_path'):
                self.DB_PATH         = paths['db_path']

        self.GENBANK_FILE    = APP_ROOT + self.GENBANK_FILE
        self.PROMOTER_FILE   = APP_ROOT + self.PROMOTER_FILE
        self.TERMINATOR_FILE = APP_ROOT + self.TERMINATOR_FILE
        self.GENOME_SEQUENCE = APP_ROOT + self.GENOME_SEQUENCE

        # self.CDS_OUT         = APP_ROOT + self.CDS_OUT
        # self.rRNA_OUT        = APP_ROOT + self.rRNA_OUT
        # self.tRNA_OUT        = APP_ROOT + self.tRNA_OUT
        # self.PROMOTER_OUT    = APP_ROOT + self.PROMOTER_OUT
        # self.TERMINATOR_OUT  = APP_ROOT + self.TERMINATOR_OUT

        self.DB_PATH         = APP_ROOT + self.DB_PATH

    def cleanup(self):
        # self.__cleanup(
        #     self.CDS_OUT, self.rRNA_OUT, self.tRNA_OUT, self.PROMOTER_OUT,
        #     self.TERMINATOR_OUT)
        pass

    def __cleanup(self, *filenames):
        for filename in filenames:
            if os.path.isfile(filename):
                os.remove(filename)

    def find_app_root(self):
        root = os.path.dirname(__file__)
        while not os.path.exists(os.path.join(root, 'setup.py')):
            root = os.path.abspath(os.path.join(root, os.path.pardir))
        return root

class Mapper(object):

    def __init__(self, db_config):
        self.conf = db_config

        if os.path.isfile(self.conf.DB_PATH):
            self.reflection = True
            # print "DB[%s] is exist" % (self.DB_PATH)
        else:
            self.reflection = False
            # print "DB[%s] is NOT exist" % (self.DB_PATH)

        self.engine = create_engine(
            'sqlite:///' + self.conf.DB_PATH, echo=False)

        if not self.reflection:
            # print "Reflection is OFF"
            species.Base.metadata.create_all(self.engine)
        elif self.reflection:
            # print "Reflection is ON"
            metadata = MetaData(bind=self.engine)
            metadata.reflect(bind=self.engine)

        maker = sessionmaker(bind=self.engine)
        self.session = maker()

    def __destroy_DB(self):
        self.session.delete()

    def generate(self):
        if not self.reflection:
            # print "Generating DB..."
            generate_decs(self.session, self.conf)
        elif self.reflection:
            # print "Use database reflection..."
            pass

class QueryBuilder(Mapper):

    def __init__(self, db_config):
        Mapper.__init__(self, db_config)

        self.generate()

    def count_stored_records(self):
        return self.session.query(species.CDSDec).filter(species.CDSDec.start).count()

    def collect_cds_records(self):
        all_rec = []
        for row in self.session.query(species.CDSDec).all():
            all_rec.append(row)
        return all_rec

    def collect_all_gene_name(self):
        names = []
        for record in self.session.query(species.CDSDec).order_by(species.CDSDec.name):
            names.append(record.name)
        return names

    def find_by_name(self, gene_name):
        for rec in self.session.query(species.CDSDec).filter_by(name=gene_name):
            return rec

    def collect_annotations_filter_by_strand(self, strand):
        filt_recs = []
        if strand == 1 or strand == -1:
            for rec in self.session.query(species.CDSDec).filter_by(strand=strand):
                filt_recs.append(rec)
            return filt_recs
        else:
            raise RuntimeError, "Invalid argument given [%d], must be -1 or 1" % (strand)

    def include_record_in_region(self, start, end):
        records = []
        for record in self.session.query(species.CDSDec).filter(species.CDSDec.start.between(start, end)):
            records.append(record)
        return records

    def include_gene_in_region(self, start, end):
        genes = []
        for gene in self.session.query(species.CDSDec).filter(species.CDSDec.start.between(start, end)):
            # genes.append(gene.name)
            genes.append(gene)
        return genes

    def include_seq_in_region(self, start, end):
        seq = []
        for s in self.session.query(species.CDSDec).filter(species.CDSDec.start.between(start, end)):
            seq.append(s.sequence)
        return seq
