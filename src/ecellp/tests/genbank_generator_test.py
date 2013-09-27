#!/usr/bin/env python

import ecellp.genbank_generator as genbank
import ecellp.session as session

import unittest


class GenbankDecGeneratorTest(unittest.TestCase):
    
    def setUp(self):
        gbk_file = ""
        self.gbk = genbank.GenbankDecGenerator(gbk_file)
        
    def tearDown(self):
        print "Done"

    def test_generate(self):
        pass

    def suite():
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(GenbankDecGeneratorTest))
        return suite

        
if __name__ == '__main__':
    unittest.main(verbosity=2).run(suite)
