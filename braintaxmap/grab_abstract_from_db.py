from py2neo import Graph

from braintaxmap.config import neo4j_db_creds, neo4j_URL
from braintaxmap.data_processing import getmnemdef

# https://stackoverflow.com/questions/37295674/how-to-get-a-dictionary-or-list-from-class-py2neo-database-record


if __name__ == '__main__':
    graph = Graph(neo4j_URL,auth=neo4j_db_creds)

    mnemonic_translate = getmnemdef()
    records = graph.run("match (n:article) return n limit 5")

    unknown_keys = []  # found that documentation of medline.record apparently is wrong?

    for rec in records:
        for k, v in rec.items()[0][1].items():
            try:
                print(mnemonic_translate[k])
                print(v)
                print(25 * '-')
            except KeyError:
                unknown_keys.append((k, v))

    input('continue to uknown keys... (press enter)')

    for k, v in unknown_keys:
        print(f'unknown: [{k}]')
        print(v)
        print(25 * '-=')
