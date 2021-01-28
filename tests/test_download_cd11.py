import os 
import re 
import unittest

from typing_extensions import final

# how to test using unittest:
# cd project
# python -m unittest tests.test_download_icd11
# use -b for supressing the prints from functions

from braintaxmap.tools import download_andunzipicd11

class TestTabulationDownloader(unittest.TestCase):
    def test_download(self):
        #setup
        DATA_DIR = os.path.join(*[os.getcwd(),'data'])
        temploc=os.path.join(*[DATA_DIR,'icd11 archive.zip'])
        finalloc=os.path.join(*[DATA_DIR,'lists-other'])
        final_filename=os.path.join(*[DATA_DIR,'lists-other','simpleTabulation.xlsx'])
    
        download_andunzipicd11(temploc=temploc, finalloc=finalloc)
        print('done download...')
        print('zipfile exists?',os.path.isfile(temploc))
        print('exlxfile exists?',os.path.isfile(final_filename))

if __name__ == '__main__':

    unittest.main()
