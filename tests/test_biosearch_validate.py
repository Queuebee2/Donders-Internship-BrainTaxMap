# check whether adding/subtracting verbs from the original list of verbs makes any difference.
# python -m unittest tests.<filename>
import unittest
from os import read

from taxmap.biosearch import *

import tests


class TestBioSearch(unittest.TestCase):
    def setUp(self):
        print("Test:", self._testMethodName)

    def test_validate_email(self):
        validate_settings()

    def test_default_email(self):
        validate_settings(email="your.devs.email@email.domain")

    def test_default_api_key(self):
        validate_settings(api_key="s0meapikey")


if __name__ == '__main__':
    unittest.main()
