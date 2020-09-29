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



"""
DATA_DIR = ".." + os.sep + 'data' + os.sep
INPUT_JSON_FILENAME = DATA_DIR+'mouse-brain-atlas_1-structure-hirarchy_1.json'
# OUTPUT_JSON_FILENAME = DATA_DIR+"unpacked_brain_atlas_reversed.pickle"



def unpack3(json):

    result = dict()

    def recursive_iter(obj, keys=()):

        if isinstance(obj, dict):
            # skip the msg node
            if 'children' in obj.keys():
                if obj['children'] == []:
                    result[obj['name']] = obj
                else:
                    copy = obj.copy()
                    copy['children'] = [c['name'] for c in copy['children']]
                    result[copy['name']] = copy

            for k, v in obj.items():
                yield from recursive_iter(v, keys + (k,))
        elif any(isinstance(obj, t) for t in (list, tuple)):
            for idx, item in enumerate(obj):
                yield from recursive_iter(item, keys + (idx,))
        else:
            yield keys, obj

    for a in recursive_iter(json):
        pass

    return result

from raymondfuzz import fuzz

def do_unpack():
    jsonobj = json.load(open(INPUT_JSON_FILENAME))
    all_nodes_in_json = unpack3(jsonobj['msg'])
    return all_nodes_in_json


if __name__ == '__main__':

    # load original mouse brain structure hirarchy json file
    jsonobj = json.load(open(INPUT_JSON_FILENAME))

    all_nodes_in_json = unpack3(jsonobj['msg'])
    c=0
    for k,v in all_nodes_in_json.items():
        print(k,v)
        # fuzz()
        if k == 'secondary motor area':
            break
        c+=1

    print(c)
