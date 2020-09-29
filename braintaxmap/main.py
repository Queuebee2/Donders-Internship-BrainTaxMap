# Author        Milain Lambers
# Github        Queuebee2
# Date          8-9-2020

""" Work in progress

From a hirarchy of structures, the end nodes are to be linked to their respective function. 

- relations will be determined/tested with literature  ( ? )


"""
import os
import sys
sys.path.append("..\\src")

from tools import raymondfuzz
from querymachine import QueryMachine

DATA_DIR = os.path.join(os.path.dirname(__file__), f'..', 'data') + os.sep
JSON_FILENAME = DATA_DIR+'mouse-brain-atlas_1-structure-hirarchy_1.json'
MESH_FILENAME = DATA_DIR+'d2020 - MESH 1.bin'


class DataNode():
    def __init__(self, id, atlas_id,ontology_id,
                acronym,name,color_hex_triplet,graph_order,
                st_level,hemisphere_id,parent_structure_id,
                children=[]):
        """ the arguments are taken from the mouse brain atlas structure """
        

class DataManager():
    def __init__(self):
        pass
    
    @staticmethod
    def readMeSH():
        """
        reads a binary file containing mesh terms
        using this as reference
        https://code.tutsplus.com/tutorials/working-with-mesh-files-in-python-linking-terms-and-numbers--cms-28587

        """
        # renames?
        with open(MESH_FILENAME, mode='rb') as file:
            mesh = file.readlines()
        
            return mesh


    @staticmethod
    def readJSON():
        """ Currently yields (unreadable!/unused) paraent nodes and a leaf node from the JSON File
        containing mouse brain structures
        """
        # rename to flatten json? readAndFlatten? NormaliseJSON?
        # should we use XML instead?
        import json
        """ Using recursion to flatten the structure?? https://towardsdatascience.com/how-to-flatten-deeply-nested-json-objects-in-non-recursive-elegant-python-55f96533103d

        not sure..."""

        jsonobj = json.load(open(JSON_FILENAME))
        
        # from here 
        # https://towardsdatascience.com/flattening-json-objects-in-python-f5343c794b10
        def flatten_json(json):
            """ Notes
                Could do with some work ....

                - check if children == [], only those are needed. (?)
            """
            out = {}

            def flatten(json, name=''):
                if type(json) is dict:
                    for a in json:
                        flatten(json[a], name + a + '_')
                elif type(json) is list:
                    i = 0
                    for a in json:
                        flatten(a, name + str(i) + '_')
                        i += 1
                else:
                    out[name[:-1]] = json

            flatten(json)
            return out
                
        data = jsonobj['msg']
        for root in data:

            # get all child-nodes per item in the 'msg',
            # this is really redundant and inefficient but the 
            # data is small so it could be acceptable? 

            children = flatten_json(root)

            # blocker for test/demo
            input('press enter to continue reading the json file')

            # the name has a weird format
            # children_0_children_1_children_0_children_1_children_2_children_6_children_3_name
            # the flattener function should be changed
            # when there's a consensus (maybe leave it, since only the node will be
            # sought for, except that we might need parent-root data later so, 
            # lots of things to consider.

            for name, v in children.items():
                if name.endswith('name'):
                    # name attribute of a end-node
                    # print(name, v)  # change to yield eventually 
                    yield name, v
                else:
                    pass



    


  


def centr_print(title='', motif='- ',amt=60):
    side = (amt//2//len(motif)) - (len(title)//2)
    print(f"{(side*motif)}{motif[0]} {title} {side*motif}")

def test_json():
    d = DataManager()
    j = d.readJSON()

    # j is flattened json generator
    # j generates/yields parentnodes(unreadable), leaf-node (a brain structure)
    for i, leaf in enumerate(j):
        print(i, leaf)

    centr_print(motif='- - - - - ')
    print(f'There are {i} leaf-nodes')
    del d,j

def test_binfile():
    d = DataManager()
    m = d.readMeSH()

    # m is the read mesh terms (in binary format still, good for regex search)
    for i in m:
        print(i)



def test_pubmed_query():
    q = QueryMachine()
        # d = DataManager()

        # m = d.readMeSH()
    m = ['mice brain', 'barrel cortex']
    """ Mouse instead mice, mesh -> singular 
    
    """

    

    for term in m:
        records = q.queryPubMed(term)
        for rec in records:
            
            rec.setdefault('MH','No MeSH terms')

            print('title:', rec['TI']  if rec['TI'] else 'No Title')
            print('MeSH Terms:', rec['MH'])
            print('Abstarct:', rec['AB']if rec['AB'] else 'No Abstract')

            # raymondfuzz.fuzz() # spread prints out ->  


def test_sqlalchemy():
    from itertools import count



    class Bob():
        _ids = count(0)

        def __init__(self):
            self.id = next(self._ids)
            self.name = f'Bob {self.id}'
            self.type = 'Bob'
            self.value = self.id*42


    bobs = [Bob() for x in range(100)]
    print(f'there are {len(bobs)} bobs')


    from sqlalchemy import Column, Integer, Unicode, UnicodeText, String
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker
    from sqlalchemy.ext.declarative import declarative_base

    from random import choice
    from string import ascii_letters as letters

    engine = create_engine('sqlite:///teste.db', echo=True)
    Base = declarative_base(bind=engine)


    class User(Base):
        __tablename__ = 'users'
        id = Column(Integer, primary_key=True)
        name = Column(Unicode(40))
        address = Column(UnicodeText, nullable=True)
        password = Column(String(20))

        def __init__(self, name, address=None, password=None):
            self.name = name
            self.address = address
            if password is None:
                password = ''.join(choice(letters) for n in xrange(10))
            self.password = password

    Base.metadata.create_all()

    Session = sessionmaker(bind=engine)
    s = Session()

def automated_search():
    print('not implemented')
    
TESTS = {"JSON":(0,test_json),
       "BINFILE":(0,test_binfile),
       "pubmedQuery":(0,test_pubmed_query),
       "SQLAlchemy":(0, test_sqlalchemy)
       }

PROGRAMS = {"AUTOMATED_ARTICLE_SEARCH":(1, automated_search)
        }

if __name__ == '__main__':

    centr_print(f'Testing stuff!')

    for testkey, (testmode, testfunc) in TESTS.items():
        if testmode:
            print(f'testing {testkey}')
            testfunc()
            centr_print(f'Done testing {testkey}')
    
    centr_print(f'All testing done')
    centr_print()
    centr_print(f'Running programs')

    
    for progkey, (mode, func) in PROGRAMS.items():
        if mode:
            print(f'running {progkey}')
            func()
            centr_print(f'Done testing {progkey}')

    centr_print(f'Done running programs')
    