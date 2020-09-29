import os
from itertools import chain, starmap
from time import sleep
from random import random
import json
import os
import pickle
import sys
"""

    should investigate this
    https://stackoverflow.com/a/38397347/6934388
    https://stackoverflow.com/questions/37315845/represent-infinite-depth-nodes-in-nested-json-form-in-neo4j

"""
DATA_DIR = ".." + os.sep + 'data' + os.sep
INPUT_JSON_FILENAME = DATA_DIR+'mouse-brain-atlas_1-structure-hirarchy_1.json'
OUTPUT_JSON_FILENAME = DATA_DIR+"extracted-json-nodes.pickle"


FUZZING = True
verbose=False

def fuzz(amt=False):
    if FUZZING:
        if not amt:
            sleep(random())
        else:
            sleep(amt)
            



def json_extract_leaf_nodes(json):
    """ Inverts/flattens a json file (specific filename in header)
    TODO rename to reverse ? invert ? 
    TODO turn into a generater function, maybe.

    """

    
    result = dict()

    for root in json['msg']:
        intermediate = dict()

        def extract_nodes(json, parent=None):

            if json['children'] == []:
                intermediate[json['name']] = json
                intermediate[json['name']]['parent'] = parent
            
            else:
                # previous parent
                p = parent
                # copy this node as the new parent for the next node
                parent = json.copy()
                # remove children from this copy
                parent['children'] = [(child['name'], child['id']) for child in json['children']]
                # set the previous parent as parent for this parent
                parent['parent'] = p
                if parent['name'] not in intermediate.keys():
                    intermediate[parent['name']] = parent
                for child in json['children']:
                    extract_nodes(child, parent)

        extract_nodes(root)
        result.update(intermediate)

    return result


def load_previous(filename=OUTPUT_JSON_FILENAME):
    """Load previously reversed json of mouse brain hirarchy file"""
    with open(filename, 'rb') as infile:
        unpacked_json = pickle.load(infile)
        print('successfully loaded')
        return unpacked_json



if __name__ == '__main__':

    # load original mouse brain structure hirarchy json file
    jsonobj = json.load(open(INPUT_JSON_FILENAME))
    sys.setrecursionlimit(10**4)
    
    try:

        # try to load previously reversed json
        print('trying to load')
        
        fuzz()
        unpacked_json = load_previous()
    except Exception as e:

        # remake a new reverse json
        print(e)
        print('an exception occurred.')
        fuzz()
        print('trying to unflatten json..')
        fuzz()
        unpacked_json = json_extract_leaf_nodes(jsonobj)

        print(len(unpacked_json))
        fuzz()

        # save the new reverse json
        with open(OUTPUT_JSON_FILENAME, 'wb') as outfile:
            pickle.dump(unpacked_json, outfile)
            print(f'succesfully saved to {OUTPUT_JSON_FILENAME}')
            fuzz(5)
    

    # print all new entries in the reversed json
    # for k, v in unpacked_json.items():
    #     print(k, v)
    #     fuzz()
    #     print()
    #     print()
    #     print()
    #     print()
    print(len(unpacked_json))
    print(unpacked_json.keys())
    

    child = unpacked_json['Orbital area, ventrolateral part, layer 1']
    for k, v in child.items():
        print('item',k, v)
                    