# check whether adding/subtracting verbs from the original list of verbs makes any difference.
# python -m unittest tests.<filename>
import unittest
from os import read

from taxmap.biosearch import *

import tests


class TestBioSearch(unittest.TestCase):
    def setUp(self):
        print("Test:", self._testMethodName)

    def test_find_ids(self):

        search_keywords = ['brain', 'amygdala']
        mesh_include = ["Rodent"]
        mesh_exclude = ["Human"]
        ids = find_ids(search_keywords,
                       mesh_include=mesh_include,
                       mesh_exclude=mesh_exclude)
        self.assertGreater(len(ids), 500000, msg='Ids cool')


if __name__ == '__main__':
    unittest.main()
