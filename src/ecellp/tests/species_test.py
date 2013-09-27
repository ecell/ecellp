#!/usr/bin/env python
# -*- coding:utf-8

import unittest

import ecellp.species as species


class SpeciesSample(unittest.TestCase):

    def test_species1(self):
        obj = species.CDSDec("test", +1, 1, 100, "cds", "ATGC")
        
        self.assertEqual(obj.name, "test")

    def test_species2(self):
        obj = species.CDSDec("test", +1, 1, 100, "cds", "ATGC")
        self.assertEqual(obj.feature, "cds")

def suite():
    suite = unittest.TestSuite()
    suite.addTests(unittest.makeSuite(SpeciesSample))
    return suite


if __name__ == '__main__':
    unittest.main(verbosity=2)
