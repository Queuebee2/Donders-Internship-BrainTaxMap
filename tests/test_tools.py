# check whether adding/subtracting verbs from the original list of verbs makes any difference.
# python -m unittest tests.<filename>
import unittest
from os import read

import tests


class TestTools(unittest.TestCase):
    def setUp(self):
        print("Test:", self._testMethodName)

    def test_backup(self):
        from taxmap.tools import backup_search
        backup_search()


if __name__ == '__main__':
    unittest.main()
