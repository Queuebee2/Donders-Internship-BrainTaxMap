import os

from querymachine import QueryMachine
from tools import centr_print, tester, programrunner

DATA_DIR = os.path.join(*[os.getcwd(),'data'])
TSV_MNEMONICS_FILENAME = os.path.join(*[DATA_DIR , 'tab_delimited_medline_record_mnemonics.tsv'])
JSON_FILENAME = os.path.join(*[DATA_DIR , 'mouse-brain-atlas_1-structure-hirarchy_1.json'])
MESH_FILENAME = os.path.join(*[DATA_DIR , 'd2020 - MESH 1.bin'])
ERROR_LOG_FILE = os.path.join(*[DATA_DIR, 'article_insert_error_log.txt'])  # todo make error output dir


def squish_tsvs():
    """ squish some data into one tabdelimited file"""
    FILE_NAMES = ["Behavioral domains and neural structures.txt",
                  "Behavioral domains and their descriptions in the literature.txt",
                  "Behavioral Domains and their neural networks.txt",
                  "Behavioral Domains, related and  task-related neural regions.txt"]

    with open(DATA_DIR + 'functional-hirarchy.tsv', 'w') as out:

        out.write("### note that after each line of ---'s, there is a row of headers\n")

        for filename in FILE_NAMES:
            title = filename[:-4]  # strips .txt
            print(title)
            out.write(f'[{title}]'.center(100, "-") + "\n")
            with open(DATA_DIR + filename, 'r') as f:
                for line in f:
                    out.write(line)


# noinspection SpellCheckingInspection
def getmnemdef():
    """ The  [Bio.Medline.Record](https://biopython.org/docs/1.75/api/Bio.Medline.html)
        class has undescriptive mnemonics instead of full string keys
        this function parses the mnemonics and their 'interpretation' into a dictionary
    """
    definitions = dict()

    with open(TSV_MNEMONICS_FILENAME, 'r') as file:
        for l in file:
            line = l.strip('\n')
            [mnemonic, definition] = line.split('\t')
            definitions[mnemonic] = definition

    return definitions


def matrixify_abstract(text, delim='-'):
    """unused : """

    useless = ['by', 'the', 'way', 'for', 'all', 'in', 'with', 'a']
    words = [w for w in text.split() if w not in useless]

    combos = dict()  # TODO consider default dict or

    for i, word in enumerate(words):
        if i > 0:  # skip first keyerror and avert taking first and lastword
            other_word = words[i - 1]
            if f'{word}{delim}{other_word}' in combos.keys():
                combos[f'{word}{delim}{other_word}'] += 1
            elif f'{other_word}{delim}{word}' in combos.keys():
                combos[f'{other_word}{delim}{word}'] += 1
            else:
                combos[f'{word}{delim}{other_word}'] = 1

    for combo, count in sorted(combos.items(), key=lambda item: combos[item[0]]):
        print(combo, count)

    return combos


def test_pubmed_query():
    q = QueryMachine()
    # d = DataManager()

    # m = d.readMeSH()
    m = ['mice brain', 'barrel cortex']
    """ Mouse instead mice, mesh -> singular 

    """

    for term in m:
        records = q.queryPubMed(term)
        for rec in records:
            rec.setdefault('MH', 'No MeSH terms')

            print('title:', rec['TI'] if rec['TI'] else 'No Title')
            print('MeSH Terms:', rec['MH'])
            print('Abstarct:', rec['AB'] if rec['AB'] else 'No Abstract')

            # raymondfuzz.fuzz() # spread prints out ->


def run_mnemdef():
    translate_mnemonics = getmnemdef()
    c = 0
    for k, v in translate_mnemonics.items():

        print(k, v)

        if c > 5:
            break
        c += 1


def automated_search():
    print('not implemented')


TESTS = {"pubmedQuery": (0, test_pubmed_query),
         'mnemonic definitions': (1, run_mnemdef),
         }

PROGRAMS = {"AUTOMATED_ARTICLE_SEARCH": (1, automated_search)
            }

if __name__ == '__main__':
    tester(TESTS)
    centr_print()
    programrunner(PROGRAMS)
