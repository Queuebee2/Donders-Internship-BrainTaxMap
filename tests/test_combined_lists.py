
# check whether adding/subtracting verbs from the original list of verbs makes any difference.
# python -m unittest tests.test_combined_lists 
from os import read
import unittest
from braintaxmap.tools import read_excluded, read_included, load_verbs

class TestTheLists(unittest.TestCase):

   
    def test_combine_lists(self):
        # not really the way to use unittests I guess
        verbs = load_verbs()
        exclude = read_excluded()
        include = read_included()
        print(f'Exclude: {len(exclude)}')
        print(f'Include: {len(include)}')
        print(f'Verbs : {len(verbs)}')
        verbs = verbs - exclude
        print(f'Verbs - excluded : {len(verbs)}')
        verbs = verbs | include
        print(f'Verbs + included: {len(verbs)}')
        

       

if __name__ == '__main__':
    unittest.main()
