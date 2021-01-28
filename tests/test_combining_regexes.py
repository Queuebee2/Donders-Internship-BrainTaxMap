
# check whether adding/subtracting verbs from the original list of verbs makes any difference.
# python -m unittest tests.test_combined_lists 
from os import read
import unittest
import re
from braintaxmap.tools import read_excluded, read_included, load_verbs, readicd11simpleTabulation
from bulkconstants import BRAIN_STRUCTURES
from sys import getsizeof
import random
class TestCombinedRegex(unittest.TestCase):

   
    def test_combine_regexes1(self):
        # not really the way to use unittests I guess
        
        disorder_hierarchy = readicd11simpleTabulation()
        disorders = set()
        for k, v in disorder_hierarchy.items():
            disorders.add(k)
            for d in v:
                disorders.add(d)
    
        verbs = load_verbs()
        
        struct_re = re.compile("|".join([rf'{word}\b' for word in BRAIN_STRUCTURES]), re.IGNORECASE)
        disorder_re = re.compile("|".join([rf'{word}\b' for word in disorders]), re.IGNORECASE)
        verbs_re = re.compile("|".join([rf'{word}\b' for word in verbs]), re.IGNORECASE)
       
        test_sentences = [
            f'{random.choice(list(BRAIN_STRUCTURES))} {random.choice(list(verbs))} {random.choice(list(disorders))}' for i in range(10)
        ]
        
        for sentence in test_sentences:
            for item in [struct_re,disorder_re,verbs_re]:
                m = item.match(sentence)
                if m:
                    print(24*'-')
                    print(sentence)
                    print('found:', m)
        
        for item in [struct_re,disorder_re,verbs_re]:
            print(getsizeof(item))

        """findingds
        this regex does not match greedily
        something short\b|something short but larger\b
        it will match 'something short' in a text with 'something short but larger' Does the order matter?
        """


       

if __name__ == '__main__':
    unittest.main()
