import os

from dotenv import find_dotenv, load_dotenv

MEDLINE_RESULTS_FILE = "Medline-results" + ".gzip"
DATA_DIR = os.path.join(os.getcwd(), 'data')
STATS_OUT_DIR = os.path.join(*[DATA_DIR, 'stats'])

load_dotenv(find_dotenv())

API_KEY = os.getenv('NCBI_API_KEY')

DISORDER_RESULT_FILENAMES = [
        #'stats_disorders found.txt',
        'all hits disorders.txt',
         #'hits all-counts.txt',
        #'stats_disorders unused.txt',
        #'hits disorders.txt',
        ]
DISORDERS = set()
VERBS = set()
BRAIN_FUNCTIONS = set()
BRAIN_STRUCTURES = set()
print('(bulkconstants.py) reading lists from files ')
for listobj, listname in [
    (DISORDERS, 'disorders'),
    (VERBS, 'verbs'),
    (BRAIN_FUNCTIONS, 'functions'),
    (BRAIN_FUNCTIONS, 'functional-hierarchy'),
    (BRAIN_STRUCTURES, 'structures')
    ]:
        items = set()
        filepath = os.path.join(DATA_DIR, 'lists', listname)
        print(f'trying..{filepath}')
        with open(filepath, 'r', encoding="utf8") as fh:
            for line in fh:
                items.add(line.strip().lower())
        listobj.update(items)
