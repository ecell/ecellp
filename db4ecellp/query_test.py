#!/usr/bin/env python 

from session import session


def query_test():

    # create Mapper class object (generating DB)
    data = '../conf.ini'
    query = session.Query(data)
    
    print query.count_stored_records()
    #print query.collect_cds_records()
    print query.find_by_name('thrL')
    #print query.collect_annotations_filter_by_strand(-1)
    print query.include_gene_in_region(21, 2012)


if __name__ == '__main__':
    query_test()
