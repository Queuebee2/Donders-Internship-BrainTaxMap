import os 
import re 
import unittest

# how to test using unittest:
# cd project
# python -m unittest tests.test_readsimpletabulation 
# use -b for supressing the prints from functions

from braintaxmap.tools import readcd11simpleTabulation

class TestTabulationReadFunction(unittest.TestCase):
    def test_readxlsxfunction(self):
        #setup
        DATA_DIR = os.getcwd() + os.sep + 'data' 
        path = os.sep.join([DATA_DIR,'lists-other','simpleTabulation.xlsx'])

        terms_found_dict = readcd11simpleTabulation(path=path)
        self.assertIsInstance(terms_found_dict, dict)

       

if __name__ == '__main__':

    unittest.main()
