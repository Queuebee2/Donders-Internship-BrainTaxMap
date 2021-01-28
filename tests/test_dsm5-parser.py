import os 
import re 
import unittest

# how to test using unittest:
# cd project
# python -m unittest tests.test_dsm5-parser -b
# use -b for supressing the prints from functions

from braintaxmap.tools import dsm5parse

class TestDSM5Parse(unittest.TestCase):
    def test_dsmparse(self):
        #setup
         DATA_DIR = os.path.join(*[os.getcwd(),'data'])
        path = os.path.join(*[DATA_DIR,'lists-other','DSM-5.txt'])

        # execute
        result = dsm5parse(path)

        self.assertIsNotNone(result)
        self.assertIsInstance(result, set)
       

if __name__ == '__main__':
    unittest.main()
