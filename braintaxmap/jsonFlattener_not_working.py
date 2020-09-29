from itertools import chain, starmap
import json
import os
import pickle

DATA_DIR = os.path.join(os.path.dirname(__file__), f'..{os.sep}..', 'data') + os.sep
JSON_FILENAME = DATA_DIR+'mouse-brain-atlas_1-structure-hirarchy_1.json'
failed_pickle = DATA_DIR+"failed_unpack_atlas.pickle"

from time import sleep
from random import random

FUZZING = True

def fuzz(amt=False):
    if FUZZING:
        if not amt:
            sleep(random())
        else:
            sleep(amt)



jsonobj = json.load(open(JSON_FILENAME))

verbose=False

def unpack(json, parents=[], final=dict()):
    if json['children'] == []: # this should be a leaf
        json['parents'] = parents
        if verbose: print('found child' + json['name'])
        return json['name'], json
    else:
        if json['name'] not in parents:
                parents.append(json['name'])
        for child in json['children']:
            if verbose: print('going deeper for' + json['name'])
            k,v = unpack(child, parents, final)
            if k == 'final':
                continue
            final[k] = v

    return 'final', final



try:
    print('trying to load')
    fuzz()
    with open(failed_pickle, 'rb') as infile:
        f = pickle.load(infile)
        print('successfully loaded')
except Exception as e:
    print(e)
    print('an exception occurred.')
    fuzz()
    print('trying to unflatten json..')
    fuzz()
    child_nodes = {}
    for item in jsonobj['msg']:
        _, f = unpack(item)
        child_nodes.update(f)
    f = child_nodes

    print(len(f))

    fuzz()
    with open(failed_pickle, 'wb') as outfile:
        pickle.dump(f, outfile)
        print('succesfully saved')
        fuzz()


c = 0
for k, v in f.items():
    print (k, v['parents'][:5], v['parents'][-5:])
    fuzz()



