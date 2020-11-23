# Author        Milain Lambers
# Github        Queuebee2
# Date          8-9-2020


""" Restructuring in progress

"""

import traceback
import logging
import pickle
import os

from querymachine import QueryMachine
from Neo4jManager import Neo4jManager

# turns a number like 2 into  2nd
ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

ERROR_LOG_FILE = 'logs' + os.sep+ 'main.py.log'

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





class searchMachine():
    id_filename = 'data' + os.sep +'found-pmids.pickle'
    search_limit = 1000000

    def __init__(self, verbose=False):

        self.logger = logging.getLogger('braintaxmap.main.searchMachine')
        
        self.db = Neo4jManager()
        self.q = QueryMachine()
        self.setup_counters()

        self.current_query = ''

        self.searches_performed = 0
        self.unique_articles_found = 0

        self.load_stored_pmids()
        self.load_MeSH_terms()
    
    def setup_counters(self):
        self.logger.info('setting all counters to 0')
        self.total_matches_found=0
        
    def load_MeSH_terms(self):
        try:
            mesh_terms = list(self.db.get_MeSH_terms())
        except Exception as e:
            self.logger.warning('Could not load mesh terms from database')
            self.logger.error(e)
            mesh_terms = set()

        finally:
            self.mesh_terms = mesh_terms

    def load_stored_pmids(self):
        """ TODO replace with get all pmids from database
        """
        try:
            with open(self.id_filename, 'rb') as infile:
                ids = pickle.load(infile)
                self.logger.info(f'loaded {len(ids)} pmids from {self.id_filename}')
        except FileNotFoundError as e:
            self.logger.warning('File for strorings found ids was was not available')
            self.logger.exception(e)
            ids = set()
        finally:
            self.stored_ids = ids

    def save_stored_pmids(self):
        with open(self.id_filename, 'wb+') as outfile:
            pickle.dump(self.stored_ids, outfile)
            self.logger.info(f'saved {len(self.stored_ids)} pmids to {self.id_filename}')


    def get_next_word(self):
        pass

    def get_unique_ids(self):
        return self.stored_ids

    def search_for_ids(self, query):

        id_list = self.q.search_pubmed_ids(query, searchlim=self.search_limit)
        self.searches_performed += 1
        self.total_matches_found += len(id_list)
        return id_list

    def filter_for_new_ids(self, id_list):

        unique_ids = set(id_list).difference(self.stored_ids)
        return unique_ids

    def log(self):
        if self.verbose:
            print(f'running query: {self.current_query}')
    
    def upload_MeSH_terms(self, new_mesh_terms):
        
        meshterms = {'name': term for term in new_mesh_terms}

        self.db.insert_new_nodes(meshterms, ['MeSH'])
    
    def scan_MeSH_terms(self, record):
        print(f"scanning mesh terms for {record['TI']}")
        self.logger.info(f"scanning mesh terms for {record['TI']}")
        try:
            mesh_terms = record['MH']

            # wonderfully analytic answer https://stackoverflow.com/a/17735466
            # TODO: consider generator expression
            # check if there are any new MeSH terms
            if not self.mesh_terms.isdisjoint(set(mesh_terms)):
                
                new_mesh_terms = set(mesh_terms).difference(self.mesh_terms)

                self.mesh_terms.update(new_mesh_terms)
                self.upload_MeSH_terms(new_mesh_terms)

            # are mesh terms relevant?
                return True
            return False

        except Exception as e:
            try:
                pub_date = record['DP']
            except Exception as err:
                pub_date = '[no publication date]'
            
            self.logger.warning(f"{record['PMID']}, from {pub_date} has no MeSH terms")
            self.logger.exception(e)
            
            return False
    

    def filter_records(self, records):
        
        for rec in records:
            
            relevant = False
            self.scan_MeSH_terms(rec)
            

    def search(self, query):
        self.current_query = query
        # get next query
        # query pubmed for ids & retrieve
        id_list = self.search_for_ids(query)

        # check which ids are new
        new_ids = self.filter_for_new_ids(id_list)
        # retrieve articles for each new i
        
        # medline records generator
        records = self.q.get_pubmed_by_pmids(list(new_ids))

        self.scan_MeSH_terms(re)

        self.filter_records(records)
    

    def run(self, runs=25):

        # threading warning (for later): this would block the program until its finished
        # but we can't do multiple searches at once!!


        for i in range(runs):
            try:
                for query in ['brain', 'neuro', 'barrel cortex']:
                    self.search(query)
                    
            except Exception as e:
                self.logger.exception(e)
            
            self.logger.error(f'main.sarchmachine.run ran again for the {ordinal(i)} time')
            


if __name__ == '__main__':

    logger.info('Starting braintaxmap main')

    s = searchMachine()

    s.run()
