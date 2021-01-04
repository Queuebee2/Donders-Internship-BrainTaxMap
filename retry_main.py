
from Bio import Entrez, Medline
from braintaxmap.config import dev_email
import nltk
print('initial email:', Entrez.email, 'initial max tries:', Entrez.max_tries)
Entrez.email = dev_email
Entrez.max_tries = 15
print('email set to', Entrez.email, 'max tries set to', Entrez.max_tries)

"""Explanation

Because sometimes, iterating over a Medline.parse(handle) generator of records, errors can occur 
that break the iteration
"""


MAX_SEARCH_LIMIT = 100000

def _testquery(query):
    """ verify query. If it's bad, raise an error"""
    if query:
        return query
    else:
        raise NotImplementedError

class PubmedSearcher():
    def __init__(self) -> None:
        self.db = 'pubmed'
        self.IDS_FILENAME='retrieved_ids.txt'
        self.ids_found=0
        pass
    
    def _clear_file(self, filename):
        open(filename, 'w').close()

    def article_loop():
        """Trying to account for all weird errors in the loop"""

        self.current_index
        while self.current_index < self.ids_found:
            pass

    def retrieve_ids(self, start=0, ):

        if start==0:
            self._clear_file(self.IDS_FILENAME)

        handle = Entrez.esearch(
            self.db, 
            term=self.query,
             email=dev_email,
             mindate='1500/01/01',
             retstart=start,
             retmax=retmax)
        result = Entrez.read(handle)
        ids = result['IdList']
        self.ids_found+=len(ids)
        if len(ids) == 0:
            return None
        print(f'found {len(ids)} ids starting with {ids[0]} ending with {ids[-1]} ')
        with open(self.IDS_FILENAME, 'a+') as out:
            for id in ids:
                out.write(f"{id}\n")
        print(f'. . . saved to file')

        return self.retrieve_ids(start+len(ids))
        

    def run(self):
        self.retrieve_ids()

    def run_query(self, query):
        self.query = _testquery(query)
        self.run()


if __name__ == '__main__':
    # fdist = nltk.FreqDist()
    # fdist[word]+=1

    p = PubmedSearcher()
    p.run_query('brain')
    print(p.ids_found)
