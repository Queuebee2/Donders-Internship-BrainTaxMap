
# Set Neo4j username and password here
# 'neo4j' is the default username.
neo4j_db_creds = ('neo4j', 'braintaxmap')

# (bolt) url to connect to db with
# to find it in Neo4j desktop go to
# your_db>manage>settings and find
# dbms.connector.bolt.listen_address=:xxxx
# put numbers xxxx in the link like so 'bolt://localhost:xxxx'
# default is 7687
neo4j_URL = 'bolt://localhost:7687'

# Used for biopython Entrez (Bio.Entrez.email)
dev_email = 'milain.lambers@gmail.com'
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
import warnings, os

DATA_DIR = ".." + os.sep + 'data' + os.sep

if neo4j_db_creds[1] == 'password':
    warnings.warn('\nWARNING\nWARNING\n\tNeo4j database username and password are still set to their defaults, '
                  'consider changing them in the braintaxmap.config.py file\nWARNING\nWARNING')

if dev_email == 'your.email@goes.here':
    warnings.warn('\nWARNING\n\\Developer email is not set. consider changing it in the braintaxmap.config.py  file\nWARNING')
