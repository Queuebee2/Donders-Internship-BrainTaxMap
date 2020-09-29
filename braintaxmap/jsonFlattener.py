from itertools import chain, starmap
from time import sleep
from random import random
from raymondfuzz import fuzz
import json
import os
import pickle
import sys

"""

    should investigate this

    https://stackoverflow.com/a/38397347/6934388



"""
DATA_DIR = ".." + os.sep + 'data' + os.sep
INPUT_JSON_FILENAME = DATA_DIR + 'mouse-brain-atlas_1-structure-hirarchy_1.json'
OUTPUT_JSON_FILENAME = DATA_DIR + "unpacked_brain_atlas_reversed.pickle"


def unpack2(json):
    def recursive_iter(obj, keys=()):
        if isinstance(obj, dict):
            for k, v in obj.items():
                yield from recursive_iter(v, keys + (k,))
        elif any(isinstance(obj, t) for t in (list, tuple)):
            for idx, item in enumerate(obj):
                yield from recursive_iter(item, keys + (idx,))
        else:
            yield keys, obj

    flatten(json)
    return result


def unpack(json):
    """ Inverts/flattens a json file (specific filename in header)
    TODO rename to reverse ? invert ?
    TODO turn into a generater function, maybe.

    """
    out = {}
    d = {'calls': 0}

    def flatten(json, parent_nodes='', leaf_node_name=False, key=False):
        if leaf_node_name and leaf_node_name != 'children':
            if leaf_node_name in out.keys():
                out[leaf_node_name][key] = json
            else:
                out[leaf_node_name] = {key: json}
                p = parent_nodes.split("_")
                while '' in p:
                    p.remove('')
                out[leaf_node_name]['parents'] = p

        else:
            if type(json) is dict:
                for key in json:
                    leaf_node_name = json['name'] if json['children'] == [] else False
                    keyname = key if key != 'children' else ''
                    flatten(json[key], parent_nodes + keyname + '', leaf_node_name, key)
                    d['calls'] += 1

            elif type(json) is list:
                i = 0
                for key in json:
                    flatten(key, parent_nodes + str(key['name']) + '_')
                    i += 1
                    d['calls'] += 1
            else:
                # parent nodes is a key at this point
                pass
                d['calls'] += 1

    flatten(json)
    print(f"unpacks in node: {d['calls']}")
    return out


def load_previous(filename=OUTPUT_JSON_FILENAME):
    """Load previously reversed json of mouse brain hirarchy file"""
    with open(filename, 'rb') as infile:
        unpacked_json = pickle.load(infile)
        print('successfully loaded')
        return unpacked_json


if __name__ == '__main__':

    # load original mouse brain structure hirarchy json file
    jsonobj = json.load(open(INPUT_JSON_FILENAME))
    sys.setrecursionlimit(10 ** 4)

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
        unpacked_json = dict()
        for root in jsonobj['msg'][0]['children']:
            print(f"unpacking {root['name']}", end=' ')
            f = unpack(root)
            unpacked_json.update(f)

        print(len(f))
        fuzz()

        # save the new reverse json
        with open(OUTPUT_JSON_FILENAME, 'wb') as outfile:
            pickle.dump(unpacked_json, outfile)
            print(f'succesfully saved to {OUTPUT_JSON_FILENAME}')
            fuzz(5)

    # print all new entries in the reversed json
    for k, v in unpacked_json.items():
        print(k)
        for s, t in v.items():
            if type(v) != list:
                print(f'{k}\t{s}\t{t}')
            else:
                for x in t:
                    print(f'{k}\t{s}\t{x}')
