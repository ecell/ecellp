#!/usr/bin/env python

import ecellp.session as session
import unittest


class QuerytoBuilderTest(unittest.TestCase):
    
    def setUp(self):
        data = '../../../conf.ini'
        self.query = session.QueryBuilder(data)
        
    def tearDown(self):
        print "Done"
        
    def test_count_stored_records(self):
        rec_count = self.query.count_stored_records()
        self.assertEqual(rec_count, 4145)

    @unittest.skip('not yet')
    def test_collect_cds_records(self):
        pass

    @unittest.skip('not yet')
    def test_collect_all_gene_name(self):
        pass

    def test_find_by_name(self):
        name = self.query.find_by_name('thrL').name
        start = self.query.find_by_name('thrL').start
        end = self.query.find_by_name('thrL').end
        strand = self.query.find_by_name('thrL').strand
        sequence = self.query.find_by_name('thrL').sequence
        feature = self.query.find_by_name('thrL').feature

        self.assertEqual(name, 'thrL')
        self.assertEqual(start, 189)
        self.assertEqual(end, 255)
        self.assertEqual(strand, 1)
        self.assertEqual(sequence,'ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA')
        self.assertEqual(feature, 'CDS')
        
    def test_collect_annotations_filter_by_strand(self):
        comp = self.query.collect_annotations_filter_by_strand(-1)
        dirc = self.query.collect_annotations_filter_by_strand(1)
        self.assertEqual(len(comp), 2130)
        self.assertEqual(len(dirc), 2015)

    @unittest.skip('not yet')
    def test_include_record_in_region(self):
        pass
        
    def test_include_gene_in_region(self):
        seq = self.query.include_gene_in_region(21, 2012)[0].values()[5]
        self.assertEqual(seq, "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA")

    @unittest.skip('not yet')
    def include_seq_ini_region(self):
        pass

    def suite():
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(QueryBuilderTest))
        return suite

        
if __name__ == '__main__':
    unittest.main(verbosity=2).run(suite)
