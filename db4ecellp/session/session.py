__program__ = 'session'
__version__ = '0.0.1'
__author__ = 'Soh Ishiguro <t10078si@sfc.keio.ac.jp>'
__copyright__ = ''
__license__ = ''

import sys
import os

from species import species
# from generators import Genbank, Promoter, Terminator, Operon, GenePromoterInteraction
from generators import Genbank, Operon, GenePromoterInteraction
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from os.path import isfile
from os import remove

import species2
import generator2


# class Mapper(Genbank, Promoter, Terminator, Operon, GenePromoterInteraction):
class Mapper(Genbank, Operon, GenePromoterInteraction):

    def __init__(self, conf):
        Genbank.__init__(self, conf)

        # if not isfile(self.CDS_OUT) and \
        #    not isfile(self.tRNA) and \
        #    not isfile(self.rRNA_OUT) and \
        #    not isfile(self.PROMOTER_OUT) and \
        #    not isfile(self.TERMINATOR_OUT):
        if True:
            self.generate_genbank_file()
            # self.generate_terminator_file()
            # self.generate_promoter_file()
            #self.generate_operon_file() #Not Implemented
            #self.generate_gene_promoter_interaction_file() #Not Implemented

        if os.path.isfile(self.DB_PATH):
            self.reflection = True
            #print "DB[%s] is exist" % (self.DB_PATH)
                        
        else:
            self.reflection = False
            #print "DB[%s] is NOT exist" % (self.DB_PATH)
        
        if not self.reflection:
            #print "Reflection is OFF"
            
            self.engine = species.create_engine('sqlite:///' + self.DB_PATH, echo=False)
            species.metadata.create_all(self.engine)
            species2.Base.metadata.create_all(self.engine)
            self.Session = species.sessionmaker()
            self.Session.configure(bind=self.engine)
            self.session = self.Session()
            
        elif self.reflection:
            #print "Reflection is ON"
            self.engine = create_engine('sqlite:///' + self.DB_PATH, echo=False)
            metadata = MetaData(bind=self.engine)
            # database reflection 
            metadata.reflect(bind=self.engine)
            self.Session = sessionmaker(bind=self.engine)
            self.session = self.Session()
            
    def __destroy_DB(self):
        self.session.delete()

    def __mapping_CDS(self):
        with open(self.CDS_OUT, 'r') as f:
            for line in f:
                (name, strand, start, end, feature, sequence) = line[:-1].split("\t")
                obj = species.CDS(name, strand, start, end, feature, sequence)
                self.session.add(obj)
        self.session.commit()

    def __mapping_tRNA(self):
        with open(self.tRNA_OUT, 'r') as f:
            for line in f:
                (name, strand, start, end, feature, sequence) = line[:-1].split("\t")
                obj = species.tRNA(name, strand, start, end, feature, sequence)
                self.session.add(obj)
        self.session.commit()

    def __mapping_rRNA(self):
        with open(self.rRNA_OUT, 'r') as f:
            for line in f:
                (name, strand, start, end, feature, sequence) = line[:-1].split("\t")
                obj = species.tRNA(name, strand, start, end, feature, sequence)
                self.session.add(obj)
        self.session.commit()

    def __mapping_promoter(self):
        gen = generator2.PromoterDecGenerator(self.PROMOTER_FILE)
        gen.generate(self.session)

    def __mapping_terminater(self):
        gen = generator2.PromoterDecGenerator(self.TERMINATOR_FILE)
        gen.generate(self.session)

    def generate_db(self):
        if not self.reflection:
            #print "Generating DB..."
            self.__mapping_CDS()
            self.__mapping_tRNA()
            self.__mapping_rRNA()
            self.__mapping_promoter()
            self.__mapping_terminater()
            
        elif self.reflection:
            pass
            #print "Use database reflection..."
            

class QueryBuilder(Mapper):

    def __init__(self, conf):
        self.conf = conf
        Mapper.__init__(self, self.conf)
        self.generate_db()
                
    def count_stored_records(self):
        return self.session.query(species.CDS).filter(species.CDS.start).count()

    def collect_cds_records(self):
        all_rec = []
        for row in self.session.query(species.CDS).all():
            all_rec.append(row)
        return all_rec
        
    def collect_all_gene_name(self):
        names = []
        for record in self.session.query(species.CDS).order_by(species.CDS.name):
            names.append(record.name)
        return names

    def find_by_name(self, gene_name):
        for rec in self.session.query(species.CDS).filter_by(name=gene_name):
            return rec
            
    def collect_annotations_filter_by_strand(self, strand):
        filt_recs = []
        if strand == 1 or strand == -1:
            for rec in self.session.query(species.CDS).filter_by(strand=strand):
                filt_recs.append(rec)
            return filt_recs
        else:
            raise RuntimeError, "Invalid argument given [%d], must be -1 or 1" % (strand)
            
    def include_record_in_region(self, start, end):
        records = []
        for record in self.session.query(species.CDS).filter(species.CDS.start.between(start, end)):
            records.append(record)
        return records

    def include_gene_in_region(self, start, end):
        genes = []
        for gene in self.session.query(species.CDS).filter(species.CDS.start.between(start, end)):
            genes.append(gene.name)
        return genes
        
    def include_seq_in_region(self, start, end):
        seq = []
        for s in self.session.query(species.CDS).filter(species.CDS.start.between(start, end)):
            seq.append(s.sequence)
        return seq

    
