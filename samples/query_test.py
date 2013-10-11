#!/usr/bin/env python

import ecellp.session
import os
import sys

def dbconf_test():
    dbconf = ecellp.session.DBConfig()
    print dbconf.DB_PATH

    if len(sys.argv) > 1:
        dbconf = ecellp.session.DBConfig(sys.argv[1])
        print sys.argv[1], dbconf.DB_PATH

def query_test():
    # create Mapper class object (generating DB)
    db_config = ecellp.session.DBConfig(os.path.abspath('./conf.ini'))
    query = ecellp.session.QueryBuilder(db_config)

    print query.count_stored_records()
    #print query.collect_cds_records()
    print query.find_by_name('thrL')
    #print query.collect_annotations_filter_by_strand(-1)
    print query.include_gene_in_region(21, 2012)


if __name__ == '__main__':
    # dbconf_test()
    query_test()
