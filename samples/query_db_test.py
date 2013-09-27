from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker

import sys
sys.path.append('../')
from species.species import Species, CDS, tRNA, rRNA, Promoter, Terminator
from session import QueryBuilder
import timeit

"""
engine = create_engine('sqlite:///test.db', echo=False)
metadata = MetaData(bind=engine)
metadata.reflect(bind=engine) # use reflection
Session = sessionmaker(bind=engine)
session = Session()
"""
#print session.query(CDS).filter_by(name="thrL").all()

def test():
    conf = '../../../conf.ini'
    q = QueryBuilder(conf)
    
    #print q.count_stored_records()
    #print q.collect_cds_records()[-1]
    #print q.collect_all_gene_name()[-1]
    print q.collect_annotations_filter_by_strand(-1)[-1]
    #print q.find_by_name("thrL")
    #print q.include_gene_in_region(50000,60000)[0]
    #print q.include_record_in_region(50000, 60000)[0]
    #print q.include_seq_in_region(2032,9059)[-1][1]

import time
t1 = time.clock()
test()
t2 = time.clock()
print t2 - t1

code = """
import sys
sys.path.append('../')
from species.species import Species, CDS, tRNA, rRNA, Promoter, Terminator
from session import QueryBuilder
conf = '../../../conf.ini'
q = QueryBuilder(conf)
print q.find_by_name("thrL")
    
"""
t = timeit.Timer(code)
print t.timeit(100)/100

