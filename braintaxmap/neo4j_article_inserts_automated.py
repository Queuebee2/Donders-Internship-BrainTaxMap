# Author    Milain Lambers
# Github    Queuebee2

import os
import pickle

import traceback

from py2neo import Graph, Node, Relationship

from config import neo4j_URL, neo4j_db_creds
from data_prep import flat_relations_hirarchy
from data_prep import load_previous
from data_processing import getmnemdef
from querymachine import QueryMachine
from tools import fuzz

DATA_DIR = ".." + os.sep + 'data' + os.sep
SOUGHT_TERM_FILENAME = DATA_DIR + 'sought_keywords.pickle'
ERROR_LOG_FILE = 'logs' + os.sep + 'article_insert_error_log.txt'  # todo make error output dir

""" see readme
this script looks for articles (not checking if they exist, yet) and inserts them into the database

"""


def load_sought(filename=SOUGHT_TERM_FILENAME, reset=False):
    """Load list of used sought_keywords"""
    if reset:
        return set()
    try:
        with open(filename, 'rb') as infile:
            sought_keywords = pickle.load(infile)
            print('successfully loaded sought_keywords')
        
        return sought_keywords
    except FileNotFoundError:
        return set()


def save_sought(sought_keywords, filename=SOUGHT_TERM_FILENAME):
    """Load list of used sought_keywords"""
    with open(filename, 'wb') as outfile:
        pickle.dump(sought_keywords, outfile)
        print('successfully dumped(saved) sought_keywords')
    
    return True


def grab_keywords():
    """ return structure  and behaviour keywords in 2 lists"""
    behaviors = iter(flat_relations_hirarchy().keys())
    structures = iter(load_previous().keys())

    return behaviors, structures

def get_custom_searchterms():
    # todo link to input field
    return [
        "((brain[Title/Abstract]) AND (behaviour[Title/Abstract])) AND (function[Title/Abstract])",
        "((brain[MeSH Major Topic]) AND (behaviour[MeSH Major Topic])) AND (function[MeSH Major Topic])",
        "barrel cortex",
        "(barrel) AND (amygdala) AND (depression)"
        "(depression) AND (hippocampus) AND (behaviour[MeSH Major Topic])"
    ]

def terms_generator(reset=False):
    """ yield keywords we haven't used yet
    
    TODO: user can insert custom keywords
    
    """

    keywords = load_sought(reset=reset)

    behaviors = iter(flat_relations_hirarchy().keys())
    structures = iter(load_previous().keys())
    custom_searchterms = iter(get_custom_searchterms)
    prev = 'none'

    for (termgenerator, labels) in [(behaviors, ['function', 'behaviour']),
                                    (structures, ['brainstructure']),
                                    (custom_searchterms, ['searched_term'])]:
        for term in termgenerator:
            while term not in keywords:
                yield term, labels
                # buffer if the program fails
                keywords.add(prev)
                save_sought(keywords) # todo  move out of this function
                prev = term

def save_MeSH_terms(record):
    pass

def harvest_articles(amt=1, reset=False, searchlimit=10000):
    """ NOTE: NO word.lower() USED YET """
    print(f'starting to harvest {amt} articles. reset_keywords_sought = {reset}')

    keywords = terms_generator(reset=reset)
    sought = load_sought()
    q = QueryMachine()

    searches_done = 0
    total_articles_found = 0
    total_unique_found = 0
    unique_pmids = set()
    unique_pmc = set()

    all_found_ids = set()
    keyword_ids_map = dict()
    pmid_pmcid_map = dict()

    for word, labels in keywords:
        fuzz(1)
        # skip already sought words
        if word in sought:
            continue


        print(f'searching for: "{word}"')
        # esearch
        id_list = set(q.search_pubmed_ids(word, searchlim=searchlimit))

        # store ids linked to this keyword
        keyword_ids_map[word] = id_list

        # keep track of total hits
        total_articles_found += len(id_list)

        # filter to get only new ids
        new_ids = set(id_list).difference(all_found_ids)

        # track total unique hits
        total_unique_found += len(new_ids)

        # efetch
        records = q.get_pubmed_by_pmids(list(new_ids))
        for rec in records:

            # todo: efficiencize
            # skips articles that are not in pmc
            # set false as default to prevent keyerror ?
            # what would be faster, keyerror > continue
            # or create False default > continue if still false
            # or check 'pmc' in keys...
            rec.setdefault('PMC', False)
            if not rec['PMC']:
                continue
            
            all_found_ids.add(rec['PMID'])
            pmid_pmcid_map[rec['PMID']] = rec['PMC']
            # set defaults for nonoccuring mnemonics
            # for mnemonic_key, interpretation in mnemonic_definitions.items():
            # rec.setdefault(mnemonic_key,'')

            # for the database, the mnemonics might not be the most readable...
            # rec = {k:v for k,v in rec.items() if v!= ''}

            unique_pmc.add(rec['PMC'])
            unique_pmids.add(rec['PMID'])

            if total_articles_found % 100 == 0:
                print(f'articles found: {total_articles_found}')
                print(f'unique pmid: {len(unique_pmids)}')
                print(f'unique pmc: {len(unique_pmc)}')
                print(60 * "-")


            yield (word, labels), rec

        searches_done += 1
        if searches_done % 2 == 0:
            print(f'searches executed: {searches_done} articles found: {total_articles_found}')

        sought.add(word)
    print(f'All searches executed, total:{searches_done}')
    print(f'Total articles found: {total_articles_found}')


def harvest_and_insert(amt=5, reset=False, searchlimit=100):
    graph = Graph(neo4j_URL, auth=neo4j_db_creds)  # more on Graph class ; https://py2neo.org/v5/database.html

    # harvester automatically gets non-used terms (if reset=False)
    # and yields new records
    article_harvester = harvest_articles(amt=amt, reset=reset, searchlimit=searchlimit)

    # definitions of mnemonic keys, e.g. 'MH' == 'Mesh Terms'
    mnemonic_definitions = getmnemdef()  # grab definitions of mnemonics

    ARTICLE_OF = Relationship.type('ARTICLE_OF')
    CITED_IN = Relationship.type('CITED_IN')  # UNUSED SO FAR.

    for (word, labels), record in article_harvester:
        word_is_parent = Node(*labels, name=word)
        word_is_parent.__primarylabel__ = labels[0]
        word_is_parent.__primarykey__ = 'name'

        article = Node('article', PMC_ID=record['PMC'], **record)
        article.__primarylabel__ = 'article'
        article.__primarykey__ = 'PMC'

        graph.merge(ARTICLE_OF(article, word_is_parent))

        print(word, record['PMC'])

        # uzz(0.05)

    print('done inserting nodes')


from datetime import datetime


def logerror(e):
    with open(ERROR_LOG_FILE, 'a') as out:
        out.write(datetime.now().strftime(
            "%I:%M%p %B %d, %Y") + " :\t" + str(e) + "\n")
        out.write(traceback.format())


if __name__ == '__main__':
    from Bio import Entrez
    from config import dev_email

    Entrez.email = dev_email

    amount = 10000
    print('starting inserts')
    runs = 0
    while runs < 25:
        try:
            harvest_and_insert(amount, reset=False, searchlimit=1000)
        except Exception as e:
            logerror(e)
            runs += 1
            fuzz()

    print('done running main')
