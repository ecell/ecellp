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

    def __init__(self, filename):
        self.__filename = filename
        if not os.path.isfile(self.__filename):
            raise RuntimeError, "Configuration file [%s] is not found" % (
                os.path.abspath(self.__filename))

        conf = ConfigParser.RawConfigParser()
        conf.read(self.__filename)
        APP_ROOT = self.find_app_root()

        # Read input file
        self.GENBANK_FILE    = APP_ROOT + conf.get('input_data', 'genbank')
        self.PROMOTER_FILE   = APP_ROOT + conf.get('input_data', 'promoter')
        self.TERMINATOR_FILE = APP_ROOT + conf.get('input_data', 'terminator')
        self.GENOME_SEQUENCE = APP_ROOT + conf.get('input_data', 'sequence')

        # Read output file
        # self.CDS_OUT        = APP_ROOT + conf.get('output_data', 'cds')
        # self.rRNA_OUT       = APP_ROOT + conf.get('output_data', 'rrna')
        # self.tRNA_OUT       = APP_ROOT + conf.get('output_data', 'trna')
        # self.PROMOTER_OUT   = APP_ROOT + conf.get('output_data', 'promoter')
        # self.TERMINATOR_OUT = APP_ROOT + conf.get('output_data', 'terminator')

        # DB
        self.DB_PATH = APP_ROOT + conf.get('db', 'db_path')

    def is_valid(self, filename):
        if not os.path.isfile(filename):
            raise IOError, "%s (wrote in %s) is not found" % (
                filename, self.__filename)

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

    def __init__(self, filename):
        self.conf = DBConfig(filename)

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

    def __init__(self, filename):
        Mapper.__init__(self, filename)

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
