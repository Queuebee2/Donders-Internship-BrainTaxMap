from braintaxmap.data_prep import load_behaviours, load_brainstructures
import os

MEDLINE_RESULTS_FILE = "Medline-results.gzip"
DATA_DIR = 'data' + os.sep + 'stats' + os.sep


# shared in multiple files

structures_pickle_path = r'data\\extracted-json-nodes.pickle'
behaviours_path = r'data\\list_of_behaviour_hirarchy.txt'
BRAIN_STRUCTURES = [n.lower()
                    for n in load_brainstructures(structures_pickle_path)]
BRAIN_FUNCTIONS = [n.lower() for n in load_behaviours(behaviours_path)]
BRAIN_STRUCTURES.remove('root')
