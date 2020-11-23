
from Bio import Entrez, Medline

from braintaxmap.config import dev_email
Entrez.email = dev_email

MAX_SEARCH_LIMIT = 100000

# todo use API key, make config file for it

class QueryMachine():
    """ 
      TODO: -set self.email in __init__
            - https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/

    ~~let's store article ID's in one file and then their stuff in another~~
    
    - mesh & fulltext from Entrez using biopython
        https://www.biostars.org/p/308345/  
    """


    def queryPMC(self, query):  # rename to get-full-text?
        """ PMC

        - E-FETCH syntax
            https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=4304705
            returns the XML of a PMC article, but we need an ID first -> look in all of PubMed first?

        - RESTFUL? 
            By default, the Europe PMC RESTful search results are sorted by relevance,
            with the most relevant result being presented first in the list. 
        
        - Module 1
            PyMed 0.8.9 https://pypi.org/project/pymed/
        """
        pass


    def querySemanticScholar(self, query):
        """ SEMANTIC SCHOLAR

        Uses AI to search through papers, categorize data, contextualise it...
        - compares citations
        - ' determine paper quality ' 
        ( source of all that  here https://www.youtube.com/watch?v=95vndf_4t4w ) 

        TODO
        

        - API KEY
            to circumvent 20/min rate limit we need an api key
            interested parties should submit a request via https://pages.semanticscholar.org/data-partners
            contact form to determine if a private API key is appropriate for your request.
        
        - GENERAL API INFO
            see more here https://api.semanticscholar.org/

        - PYTHON MODULE (?) 
            there's also a python module available https://pypi.org/project/semanticscholar/
            pip install semanticscholar

        """
        pass

   
    def queryPubMed(self, query, searchlimit=MAX_SEARCH_LIMIT, mail_adress=dev_email):  # rename to searchPubMed?
        """Search keywords against pubmed and return Bio.Medline records

        https://biopython.org/docs/1.75/api/Bio.Medline.html

        Args:
            query ([type]): [description]
        """
        if '"' not in query:
            query = f'"{query}"'
        # q.queryPubMed(term)

        # https://biopython.org/docs/1.75/api/Bio.Entrez.html#Bio.Entrez.esearch
        search_handle = Entrez.esearch(db='pubmed',
                                       sort='relevance',
                                       retmax=str(searchlimit),
                                       retmode='xml',
                                       term=query, 
                                       email=mail_adress)

        # https://biopython.org/docs/1.75/api/Bio.Entrez.Parser.html
        record = Entrez.read(search_handle)  # Bio.Entrez.Parser.DictionaryElement

        # pmids
        idlist = record["IdList"]

        fetch_handle = Entrez.efetch(db="pubmed",
                                     id=idlist,
                                     rettype="medline",
                                     retmode="text", 
                                     email=mail_adress)
        records = Medline.parse(fetch_handle)

        return records  # what should come out of here? - > ids or records

    def search_pubmed_ids(self, query, searchlim=1000000, mail_adress=dev_email):
        
        # todo: make more sophisticated way to check if a query is pre-formatted
        if '"' not in query:
            query = f'"{query}"'
        
        search_handle = Entrez.esearch(db='pubmed',
                                       sort='relevance',
                                       retmax=str(searchlim),
                                       retmode='xml',
                                       term=query, 
                                       email=mail_adress)
        
        record = Entrez.read(search_handle)  # Bio.Entrez.Parser.DictionaryElement

        return record["IdList"]
    
    def get_pubmed_by_pmids(self, pmids, mail_adress=dev_email):
        fetch_handle = Entrez.efetch(db="pubmed",
                                     id=pmids,
                                     rettype="medline",
                                     retmode="text", 
                                     email=mail_adress)

        # returns generator object with Bio.Medline.Record in it
        return Medline.parse(fetch_handle)


if __name__ == '__main__':

    # test 
    q = QueryMachine()

    ids = q.search_pubmed_ids('(brain) and (psychology) and (amygdala) and (hippocampus)')

    records = q.get_pubmed_by_pmids(ids)\

    # Bio.Medline.Record
    record = next(records)
    print('title', record['TI']) 
    print('mesh', set(record['MH']))
    print('date', record['DP'])