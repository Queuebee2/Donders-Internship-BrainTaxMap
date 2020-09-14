# Author        Milain Lambers
# Github        Queuebee2
# Date          8-9-2020

""" Work in progress

From a hirarchy of structures, the end nodes are to be linked to their respective function. 

- relations will be determined/tested with literature  ( ? )


"""

JSON_FILENAME = 'mouse-brain-atlas_1-structureGraphId_1.json'
MESH_FILENAME = 'd2020.bin'

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
        def flatten_json(y):
            """ Notes
                Could do with some work ....

                - check if children == [], only those are needed. (?)
            """
            out = {}


            def flatten(x, name=''):
                if type(x) is dict:
                    for a in x:
                        flatten(x[a], name + a + '_')
                elif type(x) is list:
                    i = 0
                    for a in x:
                        flatten(a, name + str(i) + '_')
                        i += 1
                else:
                    out[name[:-1]] = x

            flatten(y)
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



    
class QueryMachine():
    """ 

    - mesh & fulltext from Entrez using biopython
        https://www.biostars.org/p/308345/  
    """
    @staticmethod
    def queryPMC(query): # rename to get-full-text?
        """ PMC

        - E-FETCH syntax
            https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=4304705
            returns the XML of a PMC article, but we need an ID first -> look in all of PubMed first?

        - RESTFUL? 
            By default, the Europe PMC RESTful search results are sorted by relevance,
            with the most relevant result being presented first in the list. 
        
        - Module 1
            PyMed 0.8.9 https://pypi.org/project/pymed/
        """
        pass

    @staticmethod
    def querySemanticScholar(query):
        """ SEMANTIC SCHOLAR

        Uses AI to search through papers, categorize data, contextualise it...
        - compares citations
        - ' determine paper quality ' 
        ( source of all that  here https://www.youtube.com/watch?v=95vndf_4t4w ) 

        TODO
        

        - API KEY
            to circumvent 20/min rate limit we need an api key
            interested parties should submit a request via https://pages.semanticscholar.org/data-partners
            contact form to determine if a private API key is appropriate for your request.
        
        - GENERAL API INFO
            see more here https://api.semanticscholar.org/

        - PYTHON MODULE (?) 
            there's also a python module available https://pypi.org/project/semanticscholar/
            pip install semanticscholar

        """
        pass
    

    @staticmethod
    def queryPubMed(query): # rename to searchPubMed?
        """
        Search keywords against pubmed

        look at other file for example

        """

        # what should come out of here?
        pass




def test_json():
    d = DataManager()
    j = d.readJSON()

    # j is flattened json generator
    # j generates/yields parentnodes(unreadable), leaf-node (a brain structure)
    for i in j:
        print(i)

    del d,j

def test_binfile():
    d = DataManager()
    m = d.readMeSH()

    # m is the read mesh terms in (in binary format still, good for regex search)
    for i in m:
        print(i)
    


    
TESTS = {"JSON":(True,test_json),
       "BINFILE":(0,test_binfile)
       }


if __name__ == '__main__':

    print('Testing stuff')

    for testkey, (testmode, testfunc) in TESTS.items():
        if testmode:
            print(f'testing {testkey}')
            testfunc()
            print(f'Done testing {testkey}')
            print(7*"-  -  -  -  ")
    
    print('Done Testing')

    