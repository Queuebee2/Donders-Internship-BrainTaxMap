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
import textacy
import spacy

import nltk


ERROR_LOG_FILE = 'braintaxmap' + os.sep+ 'logs' + os.sep+ 'new_main.py.log'
# create logger
logger = logging.getLogger('braintaxmap.main')
logger.setLevel(logging.DEBUG)
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
from braintaxmap.tools import fuzz, load_verbs
from braintaxmap.tools import DataLoggerHelper as data

data.threshold = 1 #min amount of occurences to show

MAX_SEARCH_LIMIT = 1000000   #TODO move to new env!

nlp = spacy.load("en_core_web_lg")
CUSTOM_VERBS_SET = load_verbs('data//1000-verbs-set.txt')

def find_SVOs(sentence):
    """Try to extract SVO from a sentence and return a list of found SVO triples (s,v,o)"""
    doc = nlp(sentence)
    svo_triples = textacy.extract.subject_verb_object_triples(doc)
    svo_triples = list(svo_triples)
    return svo_triples

def find_nouns(sentence):
    """Try to find all nouns in a sentence"""
    pass

def find_verbs(sentence):
    """Try to find and extract all verbs in a sentence"""
    pass

def output_query_stats(query, stats_current, abstracts_for_output,  outfile='sample-SVO-distribution.txt'):
    """
    write stats about one query to a file
    """
    

    with open(outfile, 'a+') as fh:

        fh.write(f'-=-=-= STATS FOR {query} =-=-=-=-\n')

        data.log_topmesh(fh, query,  stats_current)
        data.log_lastmesh(fh, query,  stats_current)
        data.log_pmid_SVOs(fh, query,  stats_current)
        data.log_abstracts(fh, abstracts_for_output)
        


def output_all_stats(all_stats,  outfile='sample-SVO-distribution.txt',
    StoppedEarlyByUser=False):
    """ write stats about all queries of a run to a file
    """

    with open(outfile, 'a+') as fh:

        fh.write('\n-=-=-= TOTAL STATS -=-=-=\n')
        data.log_queries_sent(fh, all_stats)
        data.log_allmesh(fh, all_stats)
        data.log_all_svos(fh, all_stats)
        data.log_all_nounphrases(fh, all_stats)
        data.log_all_verbs(fh, all_stats)
        data.log_all_custom_verbs(fh, all_stats)

        if StoppedEarlyByUser:
            fh.write('\n !!! data incomplete, process stopped early by user !!!')


def create_datafile(outfile='sample-SVO-distribution.txt'):
    """
    ALL QUERIES
    TOP 10 ALL mesh
    TOP 10 ALL SVO

    QUERY_KEYWORD
    TOP 5 MESH
    TOP 5 SVO

    PMID
    TITLE
    ABSTRACT
    from-keywords: recorstats[pmid][found_by_keywords]
    for svo in record[SVO's]:


    SVO : [Pmids]
    
    """
    with open(outfile, 'w+') as fh:
        fh.write('top and lowest 5 abstracts per query\n\n')

def ordinal(n): 
    # turns a number like 2 into 2nd
    return "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

def scan_meshterms(record, includelist, excludelist): # or just one record?
    """ find excluded/included mesh terms in the meshterm list
     and return those lists
    """
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

def is_relevant_abstract(abstract, sets_with_relevant_terms, force_one_of_each=True):
    """within a string of text, checks whether a word from a list/set exists in it.
    can be given multiple lists. If the force_one_of_each parameter is True,
    one word from each list has to be included in the text to return True
    """
    
    checks = [0 for list in sets_with_relevant_terms]

    for word in abstract.split(' '):
        for i, s in enumerate(sets_with_relevant_terms):
            for keyword in s:
                if force_one_of_each:
                    checks[i]=1
                    break # this set is now included, skip to next
                else:
                    return True
    else:
        if force_one_of_each and all(checks):
            return True
        else:
            return False

def handle_abstract(logger, find_SVOs, DATABASE_USAGE, db, stats_all,
            stats_current, record):


    for sentence in nltk.sent_tokenize(record['AB']):
        # TODO DISCUSS LOWER() ABSTRACT?

        # handle noun phrases
        s = nlp(sentence)
        for nounphrase in s.noun_chunks:
            stats_all['nounphrases'][nounphrase.text.lower()]+=1
        
        for token in s: #TODO maybe whole ab would provide better POS on tokens
            if token.pos_ == 'VERB':
                stats_all['verbs'][str(token).lower()]+=1

        for word in nltk.word_tokenize(sentence):
            if word in CUSTOM_VERBS_SET:
                stats_all['custom_verbs'][word]+=1
        

        # handle SVO's
        svo_triples = find_SVOs(sentence)
        last_svo_id = record['PMID']
        if len(svo_triples) > 0:
            # print(svo_triples)
            pass
        for triplet in svo_triples:
            # TODO add more filters?
            if 'we' not in [str(s).lower() for s in triplet] and not (DATABASE_USAGE == 'IGNORE'):
                db.insert_svo(triplet, record['PMID'])
            logger.debug(f'added to database triplet svo: {triplet}')
        # stats_current['svo_pmids'][record['PMID']]=[str(svo) for svo in svo_triples]
        for svo in svo_triples:
            keysvo = str(svo)
            if keysvo in stats_all['allSVOs']:
                stats_all['allSVOs'][keysvo]+=1
            else:
                stats_all['allSVOs'][keysvo]=1

    return svo_triples, last_svo_id

if __name__ == '__main__':
    
    DATABASE_USAGE = 'IGNORE' # 'IGNORE', 'RETRIEVE_STORED', 'NEW_ONLY'
    
    create_datafile()
    ## inputs
    customkeywords = [] 
    queries_to_search = [
        'barrel cortex', 
        '(brain) AND (amygdala) AND (depression)',
        '(barrel cortex) AND (depression)',
    ]

    ## preparations
    ## if dotenv.get_key('NLTK_ISDOWNLOADED') else  download & dotenv.setkey....

    db = Neo4jManager()
    q = QueryMachine()

    if db.graph == None:
        logger.warning("Database is not on")
        print("Database is not on")
        exit()
    

    brainstructures = [item['node']['name'].lower() for item in db.get_brainstructures()]
    logger.info(f"retrieved {len(brainstructures)} structures")
    brainfunctions = [item['node']['name'].lower() for item in db.get_functions()]
    logger.info(f"retrieved {len(brainfunctions)} functions")
    
    queries_to_search += brainstructures[:10]

    try:
        brainstructures.remove('root ')
        brainfunctions.remove('root')
    except:
        print('no ROOT in brainfunctions?')
        pass

    relevant_terms = [set(brainstructures), set(brainfunctions)]

    # https://www.ncbi.nlm.nih.gov/mesh/68009115 Muridae
    include_MESH = ['Rats', 'Rodent', 'Mice', 'Muridae' ] # 'Brain', 'Amygdala', 'Depression'
    exclude_MESH = ['Humans' ]


    stats_all = dict(
        queries=dict(),
        allMeshtermCounts = collections.Counter(),
        allSVOs = collections.Counter(),
        svo_pmids = dict(),
        nounphrases = collections.Counter(),
        verbs=collections.Counter(),
        custom_verbs=collections.Counter()

    )
    try:
        for index_status, query in enumerate(queries_to_search):
            print(f'looking for {query}')
            
            meshterm_occurences = collections.Counter()

            filtered_pmids = dict(
                excluded_by_no_abstract = list(),
                excluded_by_irrelevant_abstract = list(),
                excluded_by_mesh = list(),
            )


            # disabled
            stats_current = dict(
                query=query,
                index=index_status,
                meshterm_occurences=meshterm_occurences,
                pmids=[],   # stats_current['pmids'] = ....
                svo_pmids=dict(),
                filtered_out=filtered_pmids,
                
            )
            
            
            query_stats = StatTracker()
            
            queryIsStructure =  query in brainstructures
            queryIsFunction = query in brainfunctions
            # search pubmed and retrieve all ids for relevant articles
            query_pmids = q.search_pubmed_ids(query, searchlim=MAX_SEARCH_LIMIT)

            stats_all['queries'][f'{query}']=query_pmids
            query_stats.add(**{f'found {len(query_pmids)} articles for query':f'{query}' })
            
            # TODO merge with StatTracker aka find better way to do stats


            # filter known pmids
            # (dont forget that at the end we still need to update the database
            #  so all found pmids link to the search-term)
            leftover_pmids = list()
            in_database_pmids = list()
            
            if DATABASE_USAGE == 'IGNORE':
                leftover_pmids = query_pmids
                records = q.get_pubmed_by_pmids(leftover_pmids)
            elif DATABASE_USAGE == 'NEW_ONLY' or 'RETRIEVE_STORED':

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

                if DATABASE_USAGE == 'RETRIEVE_STORED':
                    logger.warning('retrieving abstracts not implemented')
                    print('retrieving abstracts from db not implemented')
                    pass
                    # add records from our database to the records 
                    
            else:
                print('DATABASE_USAGE not set to valid value:')
                print('IGNORE', 'NEW_ONLY', 'RETRIEVE_STORED')
                quit()
            
            # loop over the length of leftover_pmids TODO

            # TODO See note 1 at end of file for main loop 
            # has to be adapted
            
            abstracts_for_output = []

            rec_index = 0
            while True:
                try:
                    # why try
                    """Why is everything wrapped up in a try/except?
                    """
                    record = next(records)
                    
                    rec_index+=1
            

                    if 'PMID' not in record.keys():
                        print('record without PMID? keys::')
                        for k,v in record.items():
                            print('k', k,'v', v)
                        print('Skipping record without PMID.. ?')
                        continue

                    if queryIsFunction and queryIsStructure:
                        # is there something TODO here?
                        pass
                    
                    # scan MeSHterms
                    # if only exclude-words in mesh-list, skip article
                    included_mh, excluded_mh = scan_meshterms(record, include_MESH, exclude_MESH)
                    if len(included_mh) == 0 and len(excluded_mh) > 0:
                        # skip when only excluded words exist in meshterm
                        filtered_pmids['excluded_by_mesh'].append(record['PMID'])
                        if not record['PMID'] in in_database_pmids and not (DATABASE_USAGE == 'IGNORE'): 
                            db.insert_article(record['PMID'], record, ['IRRELEVANT_MESH'])
                        continue

                    if 'MH' in record.keys():
                        # this updates the current MeSHterm counter
                        meshterm_occurences.update(record['MH'])
                        stats_all['allMeshtermCounts'].update(record['MH'])
                                                        
                    # HARD FILTERS
                    # no abstract -> out
                    if 'AB' not in record.keys():
                        filtered_pmids['excluded_by_no_abstract'].append(record['PMID'])
                        if not record['PMID'] in in_database_pmids and not (DATABASE_USAGE == 'IGNORE'): 
                            db.insert_article(record['PMID'], record, ['ABSTRACTLESS'])
                        continue

                    
                    # exclude irrelevant abstracts
                    if not is_relevant_abstract(record['AB'], relevant_terms):
                        record['isRelevant'] = '0'
                        filtered_pmids['excluded_by_irrelevant_abstract'].append(record['PMID'])
                        if not record['PMID'] in in_database_pmids and not (DATABASE_USAGE == 'IGNORE'): 
                            db.insert_article(record['PMID'], record, ['IRRELEVANT_ABSTRACT'])
                        continue

                    else:
                        
                        record['isRelevant'] = '1'
                        # TODO
                        # this is the most important part of the loop
                        svo_triples, last_svo_id = handle_abstract(logger, find_SVOs,
                         DATABASE_USAGE, db, stats_all, stats_current, record)
                                # refresh a seperate file with all stats every time
                        create_datafile('intermediate-all-stats.txt')
                        output_all_stats(stats_all, 'intermediate-all-stats.txt')
                    
                    # top 5 and lowest 5
                    topsizes=5
                    pick_this_record = ((len(leftover_pmids) - rec_index) < topsizes or \
                        (len(leftover_pmids) + rec_index) < (len(leftover_pmids) + topsizes)
                    )
            
                    if pick_this_record:
                        if record['PMID'] == last_svo_id:
                            record['SVOs'] = [str(svo) for svo in svo_triples]
                        else:
                            record['SVOs'] = ['abstract not relevant, no abstract']
                        abstracts_for_output.append(record)

                    # TODO : linked articles
                    # see note 2 about handling LINKED ARTICLES at the end of file

                    # update this record in database
                    if not record['PMID'] in in_database_pmids and not (DATABASE_USAGE == 'IGNORE'):
                        db.insert_article(record['PMID'], record, ['RELEVANT'])

                except StopIteration:
                    print(f'End of records at {rec_index} for {query}')
                    
                    break
                except Exception as e:

                    logger.exception(f'error for query "{query}"" at rec_index "{rec_index}"', e)

            query_stats.add(**{f"articles already in dtabase for f{query}":len(in_database_pmids)})
            query_stats.add(**{reason:len(value) for reason,value in filtered_pmids.items()})
            query_stats.add(**{f"MeSH term top 15 for {query}":meshterm_occurences.most_common(15)})
            query_stats.add(**{f"total # of MeSH terms for {query}":len(meshterm_occurences)})

            # stats_all['queries'].update({f'query':stats_current})

            # add stats for this query and 10 abstracts from the query top 5 and lowest 5
            output_query_stats(query, stats_current, abstracts_for_output)
            
            if not (DATABASE_USAGE == 'IGNORE'):
                for mesh_term in meshterm_occurences:
                    
                    db.insert_meshterm(mesh_term)

                        
            print(query_stats)
        
        # add total stats 
        output_all_stats(stats_all)
        


        print('done')

    except KeyboardInterrupt:
        print('program interrupted')
        print('outputting stats so far.')
        output_all_stats(stats_all)
        print('done, bye bye ♥<(^-^<)')
        quit()

    
"""Restructuring notes....


1 dealing with generators and a main loop

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

2 LINKED ARTICLES
            # if 'RF' in (record.keys()):
            #     query_stats.add(**{f"number of linked articles":len(record['RF'])})
                
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
            # abstract_tokens = nltk.word_tokenize(record['AB'])
            # tagged = nltk.pos_tag(abstract_tokens)

            # we want to strip adjectives (tagged with 'JJ.')
            # filtered_by_adjectives = [(token,tag) for token,tag in tagged if 'JJ' not in tag]
            # - (try to) stem the key words
            # - store in neo4j as is, with the verb being the relation so we get (node:<SUBJECT>)-[VERBS]-(target:<OBJECT>)
            # For example (node:brainstructure {name: "amygdala"}) - [rel:AMELIORATES] - (target:function {name: "depression"})

            # Later on, we can query the database to retrieve all different relations between nodes.
            ## TODO : short
            # add word frequency dist

            # we want to find
            # a list of meanings for each abstract for each word
            # then derive the actual relation that meaning indicates

            
            Part-of-speech (POS) Tagging	Assigning word types to tokens, like verb or noun.
            Dependency Parsing
            Entit Linking (EL)
            
            
"""


