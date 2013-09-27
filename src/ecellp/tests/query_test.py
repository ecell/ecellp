#!/usr/bin/env python

import ecellp.session as session
import unittest

class QuerytoDBTest(unittest.TestCase):
    def setUp(self):
        data = '../../../conf.ini'
        self.query = session.QueryBuilder(data)

    def tearDown(self):
        print "Done"
        
    def test_count_stored_records(self):
        
        rec_count = self.query.count_stored_records()
        self.assertEqual(rec_count, 4145)

    def test_find_by_name(self):
        gene_name = self.query.find_by_name('thrL').name
        self.assertEqual(gene_name, 'thrL')
        #print query.collect_annotations_filter_by_strand(-1)
        #print query.include_gene_in_region(21, 2012)

    def suite():
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(QuerytoDBTest))

if __name__ == '__main__':
    unittest.main(verbosity=2).run(suite)
