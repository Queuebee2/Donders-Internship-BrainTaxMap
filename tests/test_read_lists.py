import os 
import re 
import unittest

# how to test using unittest:
# cd project
# python -m unittest tests.test_dsm5-parser -b
# use -b for supressing the prints from functions

from braintaxmap.tools import readlistfile, read_included, read_excluded

class TestListReaders(unittest.TestCase):
    def test_listreader(self):
        #setup
        DATA_DIR = os.getcwd() + os.sep + 'data' 
        path = os.sep.join([DATA_DIR,'lists-to-exclude','excluded-custom-words.txt'])

        words = readlistfile(path)
        self.assertIsInstance(words, set)

    def test_excluded(self):
        words = read_excluded(verbose=False)
        self.assertIsInstance(words, set)

    def test_included(self):
        words = read_excluded(verbose=False)
        self.assertIsInstance(words, set)
       

if __name__ == '__main__':
    unittest.main()
