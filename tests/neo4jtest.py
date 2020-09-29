from py2neo import Graph, Node, Relationship
from braintaxmap.extract_nodes_from_json import load_previous
from braintaxmap.tools import fuzz

"""
https://i.stack.imgur.com/caajn.png
"""

if __name__ == '__main__':
    # connect to our graph(database)/database
    graph = Graph(password='123')
    graph.delete_all()

    tx = graph.begin()

    for label, values in load_previous().items():
        print(label, values)

        ThisNode = tx.create_unique('brainstructure', **values)

        if values['parent'] != None:
            ParentNode = tx.create_unique()


        # #anode = Node("Brainstructure", **values)
        # tx.create(anode)

        # # Then, create a second node and a relationship connecting both nodes:
        # # The zero value in the relationship tuple references the zeroth item created within that transaction, i.e. the 'parent' node.
        # parent
        # child_of = Relationship(anode, "CHILD_OF", 0)
        # parent,  = graph.create_unique(
        #     {"name":str(values['parent_structure_id'])},
        #                     child_of)

        """ maybe replace with
        create_unique(*entities)
        Create one or more unique paths or relationships in a single transaction. This is similar to create() but uses a Cypher CREATE UNIQUE clause to ensure that only relationships that do not already exist are created.
        """
        

    # tx.commit()




