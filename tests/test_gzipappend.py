import os 
import re 
import unittest
import gzip
import zipfile

# how to test using unittest:
# cd project
# python -m unittest tests.test_gzipappend 
# use -b for supressing the prints from functions

from braintaxmap.tools import dsm5parse

class TestGZIPAppend(unittest.TestCase):
    def test_gzip_appends(self):
        GZIP_TESTFILEPATH = 'gzip-test.gzip'
        with gzip.open(GZIP_TESTFILEPATH, 'ab+') as fh:
            for i in range(10):
                fh.write(f'{i}\n'.encode())

            print('peekaboo:',fh.tell())
        print('wrote')

        with gzip.open(GZIP_TESTFILEPATH, "rb") as fh:
            for i, line in enumerate(fh):
                print(i, str(line, 'utf-8'))

        
       

if __name__ == '__main__':
    unittest.main()
