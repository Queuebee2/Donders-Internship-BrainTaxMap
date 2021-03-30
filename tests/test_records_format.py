# check whether adding/subtracting verbs from the original list of verbs makes any difference.
# python -m unittest tests.test_records_format
import unittest
from pathlib import Path

from taxmap.biosearch import *

import tests


class TestRecordsFormatting(unittest.TestCase):
    def setUp(self):
        print("Test:", self._testMethodName)

    def test_records_format(self):
        from taxmap import DATA_DIR
        path_to_records = Path(DATA_DIR / 'output' / 'records' /
                               'records.txt').resolve()
        print(f"testing {path_to_records}")
        with open(path_to_records, 'r') as f:
            headers = next(f)
            self.assertEqual(headers, "PMID\tTI\tAB\tMH\n")
            for line in f:
                values = line.strip().split('\t')
                pmid = values[0]
                title = values[1]
                abstract = values[2]
                mesh = values[3]
                self.assertTrue(pmid.isnumeric())


if __name__ == '__main__':
    unittest.main()
