from py2neo import Graph

from braintaxmap.config import neo4j_db_creds, neo4j_URL

if __name__ == '__main__':
    # connect to our graph(database)/database
    graph = Graph(neo4j_URL, auth=neo4j_db_creds,)

    print('successfully connected with neo4j db')

    del graph
    
