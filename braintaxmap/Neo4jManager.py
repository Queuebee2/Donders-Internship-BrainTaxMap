
from py2neo import Graph, Node, Relationship
from braintaxmap.config import neo4j_db_creds, neo4j_URL
from json import dumps
import traceback


ARTICLE_OF = Relationship.type('ARTICLE_OF')
CITED_IN = Relationship.type('CITED_IN')  # UNUSED SO FAR.
MESHTERM_OF = Relationship.type('MESHTERM_OF') # make user decide if it connects to keyword search or only articles


import logging
db_logger = logging.getLogger('braintaxmap.Neo4jManager')


def recordsToJSON(records):
    try:
        json = dumps(records.data(), indent=4, sort_keys=True)
    except Exception as e:
        return {'error':str(e)+" stacktrace: " + traceback.format_exc()}
    return json

class ArticleNode:
    def __new__(cls, record):
        """Create py2neo.Node object and set the attributes
                the primary label is based on the PMC id. This means it is assumed the record
                describes an entry with a PMC id, (PMID is not enough)
        """
        article = Node('article', PMC_ID=record['PMC'], **record)
        article.__primarylabel__ = 'article'
        article.__primarykey__ = 'PMC'
        return article

class KeywordNode:
    def __new__(cls, word, label):
        node = Node(label, name=word)
        node.__primarylabel__ = label
        node.__primarykey__ = 'name'
        return node

class MeshNode:
    def __new__(cls, meshterm):
        node = Node('MeSH', name=meshterm)
        node.__primarylabel__ = 'MeSH'
        node.__primarykey__ = 'name'
        return node

def onNodeSelect():
    """TODO:
        - what does selecting an already expanded node do?
        - what do we do with parent/child nodes (if the node is part of a structure of its own <label> kind)
        - how do we find the node the user selected?
        - what do we need to expand out from this node?
            - do we use user-set settings to determine this
            - get the neighbouring nodes
            -

    """

    # get selected node, maybe as a param?

    # get settings to determine what to expand on

    pass


# def dummifier(func):
#     def wrapper(self, *args, **kwargs):
#         if hasattr(self, 'graph') and self.graph is not None:
#             func(self, *args, **kwargs)
#         else:
#             print(f"program tried to call '{func.__name__}({args}, {kwargs})', but there is no connection. Dummydb mode is activated")
#     return wrapper

# def for_all_methods(decorator):
#     # thanks stackoverflow
#     def decorate(cls):
#         for attr in cls.__dict__: # there's propably a better way to do this
#             if callable(getattr(cls, attr)):
#                 setattr(cls, attr, decorator(getattr(cls, attr)))
#         return cls
#     return decorate

# @for_all_methods(dummifier)

class Neo4jManager:
    def __init__(self, uri=neo4j_URL, auth=neo4j_db_creds):
        self.logger = logging.getLogger('braintaxmap.Neo4jManager.Neo4jManager')
        self.graph = None
        self.connect(uri, auth)
        
    def websiteTester(self, amount):
        cypher = f"""match (n {{name:'barrel cortex'}})-[r*0..5]-(m) return n,r,m limit {amount}"""
        result = self.run_query(cypher)
        return result

    def get_brainstructures(self):
        cypher = """Match (node:brainstructure) return node"""
        result = self.run_query(cypher)
        return result
    
    def get_behaviours(self):
        cypher = """Match (node:behaviour) return node"""
        result = self.run_query(cypher)
        return result
    
    def get_articlePMIDS(self):
        cypher = """match (node:article) return node.PMID"""
        result = self.run_query(cypher)
        return result

        
    def get_functions(self):
        cypher = """Match (node:function) return node"""
        result = self.run_query(cypher)
        return result


    def connect(self, uri, auth):
        self.logger.info('connecting to database')
        try:
            self.graph = Graph(uri, auth=auth)
            print('Neo4jManager initiated and successfully connected')
        except ConnectionRefusedError:
            self.logger.warning('connection refused')
            print('ConnectionRefusedError')
            print('either something is wrong with the authentication')
            print('or the database is not online.')
            self.graph = None
        except Exception as e:
            self.logger.error('unexpected error in Neo4jManager.Neo4jManager.connect at line 93')
            self.logger.exception(e)
            self.graph = None

    def run_query(self, query):
        try:
            results = self.graph.run(query)
        except Exception as e:
            self.logger.exception(e)
            results = []
        return results

    def insert_relation(self, child, relation_to, parent):
        """[summary]

        Args:
            child ([Create py2neo.Node]): [description]
            relation_to ([tyCreate py2neo.Relationship]): [description]
            parent ([Create py2neo.Node]): [description]
        """
        self.graph.merge(relation_to(child, parent))
        return

    def insert_svo(self, svo_triples, from_pmid):
        """ 
        insert subject verbs object into the db as nodes with 
        the verb connecting the two, attached to the article they 
        come from
        TODO; keep track of pmid in a list as attribute on node...
        at this point pmids are overwritten?
        """

        subject, verb, object = [str(s).lower() for s in svo_triples]
        verbs_relation = Relationship.type(verb.upper())


        subj_node = Node('subject', name=subject, pmid=from_pmid)
        subj_node.__primarylabel__= 'subject'
        subj_node.__primarykey__ = 'name'

        obj_node = Node('object', name=object, pmid=from_pmid)
        obj_node.__primarylabel__= 'object'
        obj_node.__primarykey__ = 'name'
        

        self.graph.merge(verbs_relation(subj_node, obj_node))

    def insert_meshterm(self, meshterm):
        """ Insert/merge meshterms into the neo4j database.

        """
        node = Node('meshterm', name=meshterm)
        node.__primarylabel__= 'meshterm'
        node.__primarykey__ = 'name'
        self.graph.merge(node)

    def insert_article(self, pmid, record, labels=['article']):
        """ Insert/merge an article into the neo4j database.
        Args:
            pmid : pubmed identifier, the primary key in the database
            record : a Bio.Medline.Record object containing article info
            labels : a list of labels for this node. Defaults to ['article']
                all labels will be changed to LOWER case and if the 'article'
                label is not present, it will be put in front of any other labels
                to serve as primary label
        """
        labels = [label.lower() for label in labels]
        labels = ['article'] + labels if 'article' not in labels else labels
        node = Node(*labels, pmid=pmid, **record)
        node.__primarylabel__= labels[0]
        node.__primarykey__ = 'PMID'
        self.graph.merge(node)
    
    def insert_node(self, nodename, attributes, labels):
        node = Node(*labels, name=nodename, **attributes)
        node.__primarylabel__ = labels[0]
        node.__primarykey__ = 'name'

        self.graph.merge(node)

    def get_all_related(self, node_name, amount=25):
        """WARNING Unsafe cypher (?)"""
        r = self.run_query(f"MATCH (start) {{ name: '{node_name}' }})--(target)"
                        f"RETURN something limit {amount}")
        return list(r)


    def get_neighbours(self, node_name, amount=25, depth=2):
        cypher = f"""
        MATCH (start {{ name: '{node_name}'}})-[r*0..{depth}]-(target) 
        RETURN start, target limit {amount}"""
        # TODO security ;/
        r = self.run_query(cypher)

        return r

    def amount_of_results(self, word):
        cypher = f"""
        MATCH (start)-[r]-(target) 
        WHERE any(prop in keys(start) where start[prop] Contains "{word}") 
        RETURN count(start) 
        LIMIT {amount}"""
        r = next(self.run_query(cypher))['count(start)']

        return r

    def find_anywhere(self, word, amount=25, depth=5):
        cypher = f"""
        MATCH (start)-[r*0..{depth}]-(target) 
        WHERE any(prop in keys(start) where start[prop] Contains "{word}") 
        RETURN start, r, target 
        LIMIT {amount}"""
        r = self.run_query(cypher)
        return r
        
    def pmid_in_database(self, pmid):
        cypher = f"""MATCH (node:article) WHERE node.PMID = "{pmid}" RETURN node"""
        r = self.run_query(cypher)
        r = list(r)
        amount = len(r)
        if amount > 1: self.logger.warning(f'more than 1 match found for pmid {pmid}')
        return amount >= 1

    def exists_in_database(self, word, label,field='name'):
        """[summary]

        Args:
            word ([str]): [the attribute value of the <field> to search for on nodes in Neo4j]
            label ([list<str>]): [list of labels, must contain at least 1]
            field (str, optional): [the field to check for the word]. Defaults to 'name'.

        Returns:
            [boolean]: [whether a node or multiple nodes exist containing <word> as a value to their
             <field> attribute]
        """        
        cypher = f"""MATCH (start:{label} {{{field}: "{word}"}}) RETURN target;"""
        r = self.run_query(cypher)
        r = list(r)
        amount = len(r)
        if amount > 1: print('more than 1 match found for', word, label, field)
        return amount >= 1
    
    def get_MeSH_terms(self):
        cypher = """Match (n:MeSH) return n"""
        result = self.run_query(cypher)
        return result
    
    def insert_new_nodes(self, node_attributes, labels):
        """insert a bunch of independent nodes with labels

        """
        for attribtues in node_attributes:
            node = Node(*labels, **attribtues)
            if 'name' in attributes:
                node.__primarykey__ == 'name'
            else:
                pass
                # TODO LOG THIS ERROR because im not sure if the 
                # set is necessary
            node.__primarylabel__ = labels[0]

            self.graph.merge(node)



    def test_random_insert(self):
        labels = ['test', 'testing', 'nostructure',]
        attributes = {'type':'testnode','value':404, 'specialPrimaryKeyTest':1337}

        node = Node(*labels, name='THIS IS A TEST NODE', **attributes)
        node.__primarylabel__ = labels[0]
        node.__primarykey__ = 'specialPrimaryKeyTest'

        self.graph.merge(node)

        cypher = 'match (n:test) delete n'
        self.run_query(cypher)

        any_tests_left = self.exists_in_database('THIS IS A TEST NODE', 'test')
        print(any_tests_left)


if __name__=='__main__':


    print('testing neo4jmanager')
    n = Neo4jManager()
    meshterms = recordsToJSON(n.get_neighbours('memory'))
    print(meshterms)

    exists = n.exists_in_database('Somatomotor areas','brainstructure')
    print('Somatomotor areas','brainstructure', 'exists?', exists)

    exists = n.exists_in_database('brain','MeSH')
    print('brain','MeSH', 'exists?', exists)

    n.test_random_insert()