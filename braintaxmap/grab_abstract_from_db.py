from py2neo import Graph
from local_neo4j_db_credentials import creds
from get_mnemonic_definition import getmnemdef

#https://stackoverflow.com/questions/37295674/how-to-get-a-dictionary-or-list-from-class-py2neo-database-record

graph = Graph(auth=creds)

mnemonic_translate = getmnemdef()
records = graph.run("match (n:article) return n limit 5")

unknown_keys = [] # found that documentation of medline.record apparently is wrong?


for rec in records:
    for k, v in rec.items()[0][1].items():
        try:
            print(mnemonic_translate[k])
            print(v)
            print(25*'-')
        except KeyError:
            unknown_keys.append((k,v))


input('continue to uknown keys... (press enter)')

for k, v in unknown_keys:
    print(f'unknown: [{k}]')
    print(v)
    print(25*'-=')