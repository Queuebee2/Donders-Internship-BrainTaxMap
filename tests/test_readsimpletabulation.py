import os 
import re 
import unittest

# how to test using unittest:
# cd project
# python -m unittest tests.test_readsimpletabulation 
# use -b for supressing the prints from functions

from braintaxmap.tools import readicd11simpleTabulation

class TestTabulationReadFunction(unittest.TestCase):
    def test_readxlsxfunction(self):
        #setup
        DATA_DIR = os.path.join(*[os.getcwd(),'data'])
        path = os.path.join(*[DATA_DIR,'lists-other','simpleTabulation.xlsx'])

        terms_found_dict = readicd11simpleTabulation(path=path)
        self.assertIsInstance(terms_found_dict, dict)

       

if __name__ == '__main__':

    unittest.main()
