from Bio import Entrez, Medline

DEVELOPER_EMAIL = "milain.lambers@gmail.com"
Entrez.email = DEVELOPER_EMAIL
# todo use API key, make config file for it

class QueryMachine():

    """ 
    let's store article ID's in one file and then their stuff in another
    
    - mesh & fulltext from Entrez using biopython
        https://www.biostars.org/p/308345/  
    """
    @staticmethod
    def queryPMC(query): # rename to get-full-text?
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

    @staticmethod
    def querySemanticScholar(query):
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
    

    @staticmethod
    def queryPubMed(query, searchlimit=10000): # rename to searchPubMed?
        """Search keywords against pubmed and return Bio.Medline records

        https://biopython.org/docs/1.75/api/Bio.Medline.html

        Args:
            query ([type]): [description]
        """
         # q.queryPubMed(term)

        # https://biopython.org/docs/1.75/api/Bio.Entrez.html#Bio.Entrez.esearch
        search_handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax=str(searchlimit),
                                retmode='xml',
                                term=query, email=DEVELOPER_EMAIL)

        # https://biopython.org/docs/1.75/api/Bio.Entrez.Parser.html
        record = Entrez.read(search_handle) # Bio.Entrez.Parser.DictionaryElement
       
        idlist = record["IdList"]

        dataList = []
        fetch_handle = Entrez.efetch(db="pubmed",
                                     id=idlist,
                                     rettype="medline",
                                     retmode="text",email=DEVELOPER_EMAIL)
        records = Medline.parse(fetch_handle)

        return records         # what should come out of here? - > ids or records
