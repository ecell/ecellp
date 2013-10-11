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

    def __init__(self, filename=None):
        '''Arguments:
            filename: OPTIONAL, indicates conf.ini path.
        '''
        self.initialize()

        if filename is not None:
            self.read(filename)

    def find_app_root(self):
        #XXX: this may cause an infinite loop
        root = os.path.dirname(__file__)
        while not os.path.exists(os.path.join(root, 'setup.py')):
            root = os.path.abspath(os.path.join(root, os.path.pardir))
        return root

    def initialize(self):
        # set default values
        self.set_paths(dict(
            GENBANK_FILE="data/NC_000913.gbk",
            PROMOTER_FILE="data/PromoterSet.txt",
            TERMINATOR_FILE="data/TerminatorSet.txt",
            GENOME_SEQUENCE="data/test.fa",
            DB_PATH="data/ecoli.sqlite3"))

    def read(self, filename):
        if not os.path.isfile(filename):
            raise (RuntimeError,
                "Configuration file [%s] is not found" % (filename))

        filename = os.path.abspath(filename)

        conf = ConfigParser.RawConfigParser()
        conf.read(filename)
        paths = {}

        if conf.has_section('input_data'):
            if conf.has_option('input_data', 'genbank'):
                paths["GENBANK_FILE"] = conf.get('input_data', 'genbank')
            if conf.has_option('input_data', 'promoter'):
                paths["PROMOTER_FILE"] = conf.get('input_data', 'promoter')
            if conf.has_option('input_data', 'terminator'):
                paths["TERMINATOR_FILE"] = conf.get('input_data', 'terminator')
            if conf.has_option('input_data', 'sequence'):
                paths["GENOME_SEQUENCE"] = conf.get('input_data', 'sequence')

        if conf.has_section('db'):
            if conf.has_option('db', 'db_path'):
                paths["DB_PATH"] = conf.get('db', 'db_path')

        self.set_paths(paths, os.path.dirname(filename))
        # self.set_paths(paths)

    def set_paths(self, paths, root_path=None):
        if root_path is None:
            # root_path = self.find_app_root() # set default
            root_path = os.getcwd()

        accepted = (
            "GENBANK_FILE", "PROMOTER_FILE", "TERMINATOR_FILE",
            "GENOME_SEQUENCE", "DB_PATH")
        for key, value in paths.items():
            if key not in accepted:
                raise RuntimeError, "unexpected key [%s] given" % (key)
            value = os.path.normpath(os.path.join(root_path, value))
            setattr(self, key, value)

    def cleanup(self):
        # self.__cleanup(
        #     self.CDS_OUT, self.rRNA_OUT, self.tRNA_OUT, self.PROMOTER_OUT,
        #     self.TERMINATOR_OUT)
        pass

    def __cleanup(self, *filenames):
        for filename in filenames:
            if os.path.isfile(filename):
                os.remove(filename)

class Mapper(object):

    def __init__(self, db_config):
        self.conf = db_config

        if os.path.isfile(self.conf.DB_PATH):
            self.reflection = True
            # print "DB[%s] is exist" % (self.DB_PATH)
        else:
            self.reflection = False
            # print "DB[%s] is NOT exist" % (self.DB_PATH)

        if not os.path.isdir(os.path.dirname(self.conf.DB_PATH)):
            os.makedirs(os.path.dirname(self.conf.DB_PATH))

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
