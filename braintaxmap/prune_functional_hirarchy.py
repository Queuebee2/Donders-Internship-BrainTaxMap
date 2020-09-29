from itertools import chain, starmap
import json
import os
import pickle

DATA_DIR = os.path.join(os.path.dirname(__file__), f'..{os.sep}..', 'data') + os.sep


"""Data from ... https://brainmap.org/taxonomy/behaviors.html

    removed all lines with -, weird A and other stuff.
    manually pruned some stuff
    removed single top hirarchy as it will be parsed anyway

    regex replaced [.]\n --> .
"""


def dot_notation_to_dict(dictionary, dots):
    prev = dictionary
    for dot in dots:
        prev.setdefault(dot, {})
        prev = prev[dot]

    return dictionary

def parse_pruned_behavioural_hirarchy():

    behaviour = dict()
    with open(IN_FILENAME, 'r') as file:

        for l in file:
            line = l.split('Â')[0] if 'Â' in l else l
            line = line.strip('\n')
            line = line.lower()


            if '.' in line:
                nodes = line.split('.')

                behaviour = dot_notation_to_dict(behaviour, nodes)

    return behaviour

def flat_relations_hirarchy():

    relations = dict()

    with open(IN_FILENAME, 'r') as file:

        for l in file:
            line = l.split('Â')[0] if 'Â' in l else l
            line = line.strip('\n')
            line = line.lower()

            if '.' in line:
                nodes = line.split('.')

                relations.setdefault(nodes[0], 'BEHAVIOR_ROOT')
                if len(nodes) == 2:
                    relations.setdefault(nodes[-1], nodes[-2])
                elif len(nodes) == 3:
                    relations.setdefault(nodes[-1], nodes[-2])
                    relations.setdefault(nodes[-2], nodes[-3])
                elif len(nodes) == 4:
                    relations.setdefault(nodes[-1], nodes[-2])
                    relations.setdefault(nodes[-2], nodes[-3])
                    relations.setdefault(nodes[-3], nodes[-4])

                else:
                    print('level 5? not accounted for')
                    raise Exception('level 5? not accounted for')



    return relations

if __name__ == "__main__":
    s = parse_pruned_behavioural_hirarchy()
    for k,v in s.items():
        if v != {}:
            for w,x in v.items():
                if x != {}:
                    for y,z in x.items():
                        print(k,w,y,z)
                else:
                    print(k,w,x)
        else:
            print(k)


