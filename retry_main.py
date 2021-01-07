
import nltk
from Bio import Entrez, Medline

from braintaxmap.config import dev_email

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
    """
    new approach: first download EVERYTHING, then analyse
    better approach ?: run downloads synchronously -> but how
    not good for web-app methinks

    p = PubmedSearcher()
    p.run_query('(some special) AND (query)')

    - stores pmids in a file -> tested with 3m results
    - iterates over pmids in file and retrieves articles in increments of <retmax>
        - should account for all errors:
            - Bio.Medline.parse or Bio.Medline.read
            - http.client.IncompleteRead(0 bytes read)
            - ConnectionResetError(10054) (forcibly closed )
            - ConnectionTimeoutError
    - when all articles are retrieved, loop over those and perform analysis

    """
    retmax = 10**7

    def __init__(self) -> None:
        self.db = 'pubmed'
        self.IDS_FILENAME = 'retrieved_ids.txt'
        self.ids_found = 0
        pass

    def _clear_file(self, filename):
        open(filename, 'w').close()

    def article_loop():
        """Trying to account for all weird errors in the loop, 
        keep exact track of articles and only download in increments of self.retmax
        downloading each article seperately would be really really
        really really slow
        maybe it'd be better to download the whole
        (gzipped) pbumed database (~130Gb)"""

        self.current_index
        while self.current_index < self.ids_found:
            pass

    def retrieve_ids(self, start=0, ):

        if start == 0:
            self._clear_file(self.IDS_FILENAME)

        handle = Entrez.esearch(
            self.db,
            term=self.query,
            email=dev_email,
            mindate='1500/01/01',
            retstart=start,
            retmax=self.retmax)

        result = Entrez.read(handle)
        ids = result['IdList']
        self.ids_found += len(ids)
        if len(ids) == 0:
            return None
        print(
            f'found {len(ids)} ids starting with {ids[0]} ending with {ids[-1]} ')
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
