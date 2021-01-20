import os 
import re 
import unittest

from typing_extensions import final

# how to test using unittest:
# cd project
# python -m unittest tests.test_download_cd11
# use -b for supressing the prints from functions

from braintaxmap.tools import download_andunzipcd11

class TestTabulationDownloader(unittest.TestCase):
    def test_download(self):
        #setup
        DATA_DIR = os.getcwd() + os.sep + 'data' 
        temploc=os.sep.join([DATA_DIR,'cd11 archive.zip'])
        finalloc=os.sep.join([DATA_DIR,'lists-other'])
        final_filename=os.sep.join([DATA_DIR,'lists-other','simpleTabulation.xlsx'])
    
        download_andunzipcd11(temploc=temploc, finalloc=finalloc)
        print('done download...')
        print('zipfile exists?',os.path.isfile(temploc))
        print('exlxfile exists?',os.path.isfile(final_filename))

if __name__ == '__main__':

    unittest.main()
