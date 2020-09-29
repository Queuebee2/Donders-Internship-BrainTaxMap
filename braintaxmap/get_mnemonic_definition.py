from itertools import chain, starmap
import json
import os
import pickle

""" The  [Bio.Medline.Record](https://biopython.org/docs/1.75/api/Bio.Medline.html)
    class has undescriptive mnemonics instead of full string keys
    this script parses the mnemonics and their 'interpretation' into a dictionary
"""
DATA_DIR = os.path.join(os.path.dirname(__file__), f'..{os.sep}..', 'data') + os.sep
TSV_MNEMONICS_FILENAME = DATA_DIR+'tab_delimited_medline_record_mnemonics.tsv'


def getmnemdef():

    definitions = dict()

    with open(TSV_MNEMONICS_FILENAME, 'r') as file:
        for l in file:
            line = l.strip('\n')
            [mnemonic, definition] = line.split('\t')
            definitions[mnemonic] = definition

    return definitions


if __name__ == '__main__':
    
    menmonic_definitions =  getmnemdef()
    for k,v in menmonic_definitions.items():
        print(f"key[{k}] has item [{v}]")