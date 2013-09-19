__program__ = 'session'
__version__ = '0.0.1'
__author__ = 'Soh Ishiguro <t10078si@sfc.keio.ac.jp>'
__copyright__ = ''
__license__ = ''

from species import species
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
import os.path
from os import remove

class Mapper(object):
    def __init__(self):
        db_path = '../../db/ecoli.sqlite3'
        if os.path.isfile(db_path):
            os.remove(db_path)
            
        self.engine = species.create_engine('sqlite:///' + db_path, echo=False)
        species.metadata.create_all(self.engine)
        self.Session = species.sessionmaker()
        self.Session.configure(bind=self.engine)
        self.session = self.Session()

    def __destroy_DB(self):
        self.session.delete()

    def __mapping_CDS(self):
        with open('../../data/CDS_annotation.tbl', 'r') as f:
            for line in f:
                (name, strand, start, end, feature, sequence) = line[:-1].split("\t")
                obj = species.CDS(name, strand, start, end, feature, sequence)
                self.session.add(obj)
        self.session.commit()

    def __mapping_tRNA(self):
        with open('../../data/tRNA_annotation.tbl', 'r') as f:
            for line in f:
                (name, strand, start, end, feature, sequence) = line[:-1].split("\t")
                obj = species.tRNA(name, strand, start, end, feature, sequence)
                self.session.add(obj)
        self.session.commit()

    def __mapping_rRNA(self):
        pass

    def __mapping_promoter(self):
        pass

    def __mapping_terminater(self):
        pass
    
    def generate_db(self):
        self.__mapping_CDS()
        self.__mapping_tRNA()


class Query(Mapper):

    def __init__(self):
        Mapper.__init__(self)
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

    
