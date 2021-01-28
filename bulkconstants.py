import os

from dotenv import find_dotenv, load_dotenv

from braintaxmap.data_prep import load_behaviours, load_brainstructures

MEDLINE_RESULTS_FILE = "Medline-results" + ".gzip"
DATA_DIR = os.path.join(os.getcwd(), 'data')
STATS_OUT_DIR = os.path.join(*[DATA_DIR, 'stats'])


load_dotenv(find_dotenv())

API_KEY = os.getenv('NCBI_API_KEY')


structures_pickle_path = os.path.join(
    *[DATA_DIR, 'extracted-json-nodes.pickle'])
behaviours_path = os.path.join(*[DATA_DIR, 'list_of_behaviour_hirarchy.txt'])
BRAIN_STRUCTURES = [n.lower()
                    for n in load_brainstructures(structures_pickle_path)]
BRAIN_FUNCTIONS = [n.lower() for n in load_behaviours(behaviours_path)]
BRAIN_STRUCTURES.remove('root')

# todo : move other lists here aswell, maybe create a bulkprep with all the prep/loading functions for these
# - icd11, dsm5, verbs..
