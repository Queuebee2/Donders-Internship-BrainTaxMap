# Author    Milain Lambers
# Github    Queuebee2

from py2neo import Graph, Node, Relationship
from prune_functional_hirarchy import flat_relations_hirarchy
from raymondfuzz import fuzz



def insertNodes():
    graph = Graph(password='123')   # more on Graph class ; https://py2neo.org/v5/database.html
    # # # # graph.# delete_all()    # DANGEROUS. ONLY USE WHEN TESTING. DELETES EVERYTHING.


    # transaction = graph.begin()

    behavioural_hirarchy = flat_relations_hirarchy() # loads previous json
    # graph.schema.create_uniqueness_constraint('function', 'name')
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


if __name__ == '__main__':

    # insertNodes()
    print('already inserted nodes')



