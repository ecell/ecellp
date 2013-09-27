#!/usr/bin/env python

import ecellp.session


def query_test():
    # create Mapper class object (generating DB)
    data = './conf.ini'
    query = ecellp.session.QueryBuilder(data)

    print query.count_stored_records()
    #print query.collect_cds_records()
    print query.find_by_name('thrL')
    #print query.collect_annotations_filter_by_strand(-1)
    print query.include_gene_in_region(21, 2012)


if __name__ == '__main__':
    query_test()
