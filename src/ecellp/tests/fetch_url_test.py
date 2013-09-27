import utils
import unittest

class FetchUrlTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_fetch_url(self):
        file = "test.out"
        url = ("http://regulondb.ccg.unam.mx/"
            "menu/download/datasets/files/PromoterSet.txt")
        result = utils.fetch_url(url, file)
        self.assertEqual(result, False)

    def suite():
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(FetchUrlTest))
        return suite


if __name__ == '__main__':
    unittest.main(verbosity=2).run(suite)
