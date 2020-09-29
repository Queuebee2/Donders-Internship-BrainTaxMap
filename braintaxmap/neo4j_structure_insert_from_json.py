# Author    Milain Lambers
# Github    Queuebee2

from py2neo import Graph, Node, Relationship
from extract_nodes_from_json import load_previous
from raymondfuzz import fuzz
from recursive_iter_test import unpack3


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


def insertNodes():
    graph = Graph(password='123')   # more on Graph class ; https://py2neo.org/v5/database.html
    graph.delete_all()              # DANGEROUS. ONLY USE WHEN TESTING. DELETES EVERYTHING.


    # transaction = graph.begin()

    structural_hirarchy = load_previous() # loads previous json
    # graph.schema.create_uniqueness_constraint('brainstructure', 'name')
    CHILD_OF = Relationship.type('CHILD_OF')

    for node_name, node_attributes in structural_hirarchy.items():

        children = node_attributes['children']
        try:
            parent_name = node_attributes['parent']['name']
        except TypeError:
            parent_name = 'Root'
        attributes = node_attributes

        del attributes['children'] # children will be put into db relations
        del attributes['parent'] # parent will be put in through relation

        parent = Node('brainstructure', name=parent_name)
        parent.__primarylabel__ = 'brainstructure'
        parent.__primarykey__ = 'name'
        this_node = Node('brainstructure', **attributes)
        this_node.__primarylabel__ = 'brainstructure'
        this_node.__primarykey__ = 'name'

        print('parent:',parent)
        print('this',this_node)

        graph.merge(CHILD_OF(this_node, parent))


        print(node_name)
        print(children)

        # """ this will be very redundant, but for 1300 entries it should be fine"""

        # print(parent)
        # # parent.__primarylabel = 'brainstructure'
        # # parent.__primarykey__ = "name"
        # # graph.merge(parent)
        # if children:
        #     for child_attributes in children:
        #         name = child_attributes[0] # id = [1]
        #         # 'lets try with name only, see if attributes get added in the end when
        #         # the child has no children'
        #         c = Node('brainstructure', name=name, )
        #         # c.__primarylabel__ = 'brainstructure'
        #         # c.__primarykey__ == 'name'
        #         # graph.merge(c)

        #         fuzz()


    print('done inserting nodes')


if __name__ == '__main__':

    # insertNodes()
    print('already inserted nodes')
    # connect to our graph(database)/database



     # transaction.commit()

    # for label, values in load_previous().items():
    #     print(label, values)

    #     ThisNode = tx.create_unique('brainstructure', **values)

    #     if values['parent'] != None:
    #         ParentNode = tx.create_unique()


    #     # #anode = Node("Brainstructure", **values)
    #     # tx.create(anode)

    #     # # Then, create a second node and a relationship connecting both nodes:
    #     # # The zero value in the relationship tuple references the zeroth item created within that transaction, i.e. the 'parent' node.
    #     # parent
    #     # child_of = Relationship(anode, "CHILD_OF", 0)
    #     # parent,  = graph.create_unique(
    #     #     {"name":str(values['parent_structure_id'])},
    #     #                     child_of)

""" maybe replace with
    #     create_unique(*entities)
        Create one or more unique paths or relationships in a single transaction. This is similar to create() but uses a Cypher CREATE UNIQUE clause to ensure that only relationships that do not already exist are created.
"""


    #




