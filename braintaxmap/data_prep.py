# Github    Queuebee2
# Author    Milain Lambers

import json
import os
import pickle
import sys

from py2neo import Graph, Node, Relationship

from braintaxmap.config import neo4j_db_creds, neo4j_URL
from braintaxmap.data_prep import load_previous
from braintaxmap.tools import centr_print, tester, programrunner, fuzz

DATA_DIR = ".." + os.sep + 'data' + os.sep
IN_FILENAME = DATA_DIR + 'list_of_behaviour_hirarchy.txt'
INPUT_JSON_FILENAME = DATA_DIR + 'mouse-brain-atlas_1-structure-hirarchy_1.json'
OUTPUT_JSON_FILENAME = DATA_DIR + "extracted-json-nodes.pickle"


class DataPrepper():
    """
    Prepares data from inputfiles and starts inserting the first
    nodes into a neo4j db

    TODO notes
        - something to save state (remember whether nodes were inserted or check db)
    """

    def __init__(self):
        name = 'dataprepper'

    def run(self, reset_db=False, delete_all=False):
        self.insert_structural_nodes(reset_db=reset_db, delete_all=delete_all)
        self.insert_behavioural_nodes()

    @classmethod
    def insert_structural_nodes(reset_db=False, delete_all=False):
        """ work in progress

        NEO4J Browser commands
            Find all nodes:
                MATCH (n:brainstructure) RETURN n
            See Scheme
                :scheme
            See constraints
                CALL db.constraints

        PURPOSE
            Insert structural (and later functional) hirarchy from json file into
            neo4j database using py2neo

        RESOURCES
            https://medium.com/neo4j/5-tips-tricks-for-fast-batched-updates-of-graph-structures-with-neo4j-and-cypher-73c7f693c8cc

            results = graph.run("MATCH (a:Person) RETURN a.name, a.born LIMIT 4").data()

        """

        # todo here check if we didn't already insert those

        graph = Graph(neo4j_URL, auth=neo4j_db_creds)  # more on Graph class ; https://py2neo.org/v5/database.html

        # starting over?
        if reset_db and delete_all:
            graph.delete_all()

        # graph.schema.create_uniqueness_constraint('brainstructure', 'name') # todo test this constraint
        structural_hirarchy = load_previous()

        CHILD_OF = Relationship.type('CHILD_OF')

        for node_name, node_attributes in structural_hirarchy.items():

            children = node_attributes['children']
            try:
                parent_name = node_attributes['parent']['name']
            except TypeError:
                parent_name = 'Root'

            attributes = node_attributes

            del attributes['children']  # children will be put into db relations
            del attributes['parent']  # parent will be put in through relation

            parent = Node('brainstructure', name=parent_name)
            parent.__primarylabel__ = 'brainstructure'
            parent.__primarykey__ = 'name'
            this_node = Node('brainstructure', **attributes)
            this_node.__primarylabel__ = 'brainstructure'
            this_node.__primarykey__ = 'name'

            # merge creates a node if it does not already exist.
            graph.merge(CHILD_OF(this_node, parent))

            print(node_name)
            print(children)


        print('done inserting nodes')

    @classmethod
    def insert_behavioural_nodes(self):
        """ see readme
        this script looks for articles (not checking if they exist, yet) and inserts them into the database

        """

        # todo here check if we didn't already insert those

        graph = Graph(neo4j_URL, auth=neo4j_db_creds)  # more on Graph class ; https://py2neo.org/v5/database.html


        behavioural_hirarchy = flat_relations_hirarchy()  # loads previous json

        # todo check on: graph.schema.create_uniqueness_constraint('function', 'name')
        SUBFUNCTION_OF = Relationship.type('SUBFUNCTION_OF')

        for function, parent in behavioural_hirarchy.items():
            parent_function = Node(*['function', 'behaviour'], name=parent)
            parent_function.__primarylabel__ = 'function'
            parent_function.__primarykey__ = 'name'

            function = Node(*['function', 'behaviour'], name=function)
            function.__primarylabel__ = 'function'
            function.__primarykey__ = 'name'

            graph.merge(SUBFUNCTION_OF(function, parent_function))

            print(function, parent)

        print('done inserting nodes')


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
    """Data from ... https://brainmap.org/taxonomy/behaviors.html

        removed all lines with -, weird A and other stuff.
        manually pruned some stuff
        removed single top hirarchy as it will be parsed anyway

        regex replaced [.]\n --> .
    """
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


def test_json_extractor():
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
        print('item', k, v)


def test_pruner():
    s = parse_pruned_behavioural_hirarchy()
    for k, v in s.items():
        if v != {}:
            for w, x in v.items():
                if x != {}:
                    for y, z in x.items():
                        print(k, w, y, z)
                else:
                    print(k, w, x)
        else:
            print(k)


TESTS = {"pruning behavioural hirarchy": (0, test_pruner),
         "extract leafnodes from json": (0, test_json_extractor),
         }
PROGRAMS = {}
if __name__ == '__main__':
    tester(TESTS)
    centr_print()
    programrunner(PROGRAMS)
