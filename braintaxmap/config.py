import os
from dotenv import load_dotenv, find_dotenv


load_dotenv(find_dotenv())

neo4j_db_id = os.getenv('NEO4J_DB_ID')
neo4j_db_pw = os.getenv('NEO4J_DB_PW')
neo4j_db_creds = (neo4j_db_id, neo4j_db_pw)
neo4j_URL = os.getenv('NEO4J_URL')
dev_email = os.getenv("DEV_EMAIL")


import warnings, os

DATA_DIR = ".." + os.sep + 'data' + os.sep

if neo4j_db_creds[1] == 'password':
    warnings.warn('\nWARNING\nWARNING\n\tNeo4j database username and password are still set to their defaults, '
                  'consider changing them in the braintaxmap.config.py file\nWARNING\nWARNING')

if dev_email == 'your.email@goes.here':
    warnings.warn('\nWARNING\n\\Developer email is not set. consider changing it in the braintaxmap.config.py  file\nWARNING')
