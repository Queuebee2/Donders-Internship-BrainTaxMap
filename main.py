# Author        Milain Lambers
# Github        Queuebee2
# Date          8-9-2020

""" ongrowing TODO list
https://stackoverflow.com/questions/1092531/event-system-in-python

more noteworthy
 The main difference is that spaCy is integrated and opinionated. 
 spaCy tries to avoid asking the user to choose between multiple algorithms that deliver equivalent functionality.
"""
import collections
import traceback
import logging
import pickle
import os
from urllib.error import HTTPError
import dotenv
from braintaxmap.config import dev_email
from urllib.error import HTTPError

import nltk


ERROR_LOG_FILE = 'braintaxmap' + os.sep+ 'logs' + os.sep+ 'new_main.py.log'
# create logger
logger = logging.getLogger('braintaxmap.main')
logger.setLevel(logging.INFO)
# create file handler which logs even debug messages
fh = logging.FileHandler(ERROR_LOG_FILE)
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)


from braintaxmap.querymachine import QueryMachine
from braintaxmap.Neo4jManager import Neo4jManager
from braintaxmap.stat_tracker import StatTracker
from braintaxmap.tools import fuzz

MAX_SEARCH_LIMIT = 100000   #TODO move to new env!





def ordinal(n): 
    # turns a number like 2 into 2nd
    return "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

def scan_meshterms(record, includelist, excludelist): # or just one record?
    included, excluded = list(), list()
    if 'MH' in record.keys():
        for m in record['MH']:
            meshterm = m.lower()
            for incl in includelist:
                if incl.lower() in meshterm:
                    included.append(meshterm)
            
            for excl in excludelist:
                if excl.lower() in meshterm:
                    excluded.append(meshterm)

    return included, excluded

def is_relevant_abstract(abstract, keywords):
    """checks if any words from the abstracts are contained in the set of keywords"""
    for word in abstract.split(' '):
        if word in keywords:
            return True
    else:
        return False

if __name__ == '__main__':


    ## inputs
    customkeywords = [] 
    queries_to_search = [
        'barrel cortex', 
        '(brain) AND (amygdala) AND (depression)',
        '(barrel cortex) AND (depression)',
        '(heartbreak)'
    ]



    ## preparations
    ## if dotenv.get_key('NLTK_ISDOWNLOADED') else  download & dotenv.setkey....

    db = Neo4jManager()
    q = QueryMachine()
    if db.graph == None:
        logger.warning("Database is not on")
        print("Database is not on")
        exit()
    

    brainstructures = [item['node']['name'] for item in db.get_brainstructures()]
    logger.info(f"retrieved {len(brainstructures)} structures")
    brainfunctions = [item['node']['name'] for item in db.get_functions()]
    logger.info(f"retrieved {len(brainfunctions)} functions")

    relevant_terms = set(brainstructures+brainfunctions)

    # https://www.ncbi.nlm.nih.gov/mesh/68009115 Muridae
    include_MESH = ['Rats', 'Rodent', 'Mice', 'Muridae' ] # 'Brain', 'Amygdala', 'Depression'
    exclude_MESH = ['Humans' ]



    for index_status, query in enumerate(queries_to_search):
        queryIsStructure =  query in brainstructures
        queryIsFunction = query in brainfunctions
        meshterm_occurences = collections.Counter()

        query_stats = StatTracker()
        # search pubmed and retrieve all ids for relevant articles
        query_pmids = q.search_pubmed_ids(query)

        query_stats.add(**{f'found {len(query_pmids)} articles for query':f'{query}' })
        
        # TODO merge with StatTracker aka find better way to do stats
        filtered_pmids = dict(
            excluded_by_no_abstract = list(),
            excluded_by_irrelevant_abstract = list(),
            excluded_by_mesh = list(),
        )

        # filter known pmids
        # (dont forget that at the end we still need to update the database
        #  so all found pmids link to the search-term)
        leftover_pmids = list()
        in_database_pmids = list()
        for pmid in query_pmids:
            if not db.pmid_in_database(pmid):
                leftover_pmids.append(pmid)
            else:
                in_database_pmids.append(pmid)

        if len(leftover_pmids) == 0:
            # no new pmids, skip query -- ? or retreive from database..? 
            logger.info(f'no new articles for {query}' )
            continue



        records = q.get_pubmed_by_pmids(leftover_pmids)
        # TODO: records should be a queue of items and the items should be removed
        # when done....
        # problem: records is a generator object, can't append.
        # un-generatoring it could mean thousands of articles in memory
        # this could be replaced by using a loop that calls next() on the generator and adding items from it
        # to a deque/queue. this way we could add items from within the loop as well
        # finding = True
        # counter=0
        # while finding:
        #     record = None
        #     try:
        #         record = next(records):
        #     except StopIteration:
        #         pass
        #     if not record:
        #         try:
        #             record = . . . next(record from linked records collection generator:?)
        #         except StopIteration:
        #             finding= False
        #      do stuff
        # annoying.
        for record in records:
            
            if queryIsFunction and queryIsStructure:
                # is there something TODO here?
                pass
            
            # scan MeSHterms
            # if only exclude-words in mesh-list, skip article
            included_mh, excluded_mh = scan_meshterms(record, include_MESH, exclude_MESH)
            if len(included_mh) == 0 and len(excluded_mh) > 0:
                # skip when only excluded words exist in meshterm
                filtered_pmids['excluded_by_mesh'].append(record['PMID'])
                if not record['PMID'] in in_database_pmids: 
                    db.insert_article(record['PMID'], record, ['IRRELEVANT_MESH'])
                continue

            if 'MH' in record.keys():
                # this updates the current MeSHterm counter
                meshterm_occurences.update(record['MH'])
                                                
            # HARD FILTERS
            # no abstract -> out
            if 'AB' not in record.keys():
                filtered_pmids['excluded_by_no_abstract'].append(record['PMID'])
                if not record['PMID'] in in_database_pmids: 
                    db.insert_article(record['PMID'], record, ['ABSTRACTLESS'])
                continue
            
            # exclude irrelevant abstracts
            if not is_relevant_abstract(record['AB'], relevant_terms):
                filtered_pmids['excluded_by_irrelevant_abstract'].append(record['PMID'])
                if not record['PMID'] in in_database_pmids: 
                    db.insert_article(record['PMID'], record, ['IRRELEVANT_ABSTRACT'])
                continue
            
            if 'RF' in (record.keys()):
                query_stats.add(**{f"number of linked articles":len(record['RF'])})
                
            # LINKED ARTICLES
            ## If the article is relevant... continue searching their linked articles?
            # that is a form of recursion... how to implement?
            # grab the record objects of each PMID we dont have yet and add it to the queue
            # (removes recursion aspect and ensures we don't get stuck in a loop)
            
            # linked_articles_pmids = q.get_pubmed_linked_pmids(record['PMID'])
            # additional_articles_to_retrieve=list()
            # for pmid in linked_articles_pmids:
            #     if not db.pmid_in_database(pmid): # AND NOT IN ADDITIONALLY_RETRIEVED GLobally # TODO
            #         additional_articles_to_retrieve.append(pmid)
            
            # linked_records = q.get_pubmed_by_pmids(additional_articles_to_retrieve)
            # for r in linked_records:
            #     records.append(r)
            # TODO: AttributeError 'generator' object has no attribute 'append'

            # https://pubmed.ncbi.nlm.nih.gov/24811994/ part of speech tagger for biomedical application
            abstract_tokens = nltk.word_tokenize(record['AB'])
            tagged = nltk.pos_tag(abstract_tokens)

            # we want to strip adjectives (tagged with 'JJ.')
            filtered_by_adjectives = [(token,tag) for token,tag in tagged if 'JJ' not in tag]

            ## TODO : short
            # add word frequency dist

            # we want to find
            # a list of meanings for each abstract for each word
            # then derive the actual relation that meaning indicates

            """
            Part-of-speech (POS) Tagging	Assigning word types to tokens, like verb or noun.
            Dependency Parsing
            Entity Linking (EL)
            """
            #
            # update this record in database
            if not record['PMID'] in in_database_pmids:
                db.insert_article(record['PMID'], record, ['RELEVANT'])    

        query_stats.add(**{"articles already in dtabase for {query}":len(in_database_pmids)})
        query_stats.add(**{reason:len(value) for reason,value in filtered_pmids.items()})
        query_stats.add(**{f"MeSH term top 10 for {query}":meshterm_occurences.most_common(10)})
        query_stats.add(**{f"total # of MeSH terms for {query}":len(meshterm_occurences)})
                
        print(query_stats)
    
    print('done')

    



