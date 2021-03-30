import gzip
import logging
import time
from collections import Counter
from pathlib import Path
from typing import Set
from urllib.error import HTTPError

import Bio.Entrez.Parser
from Bio import Entrez, Medline

from ..tools import save_list
from . import DEFAULT_RECORDS_DIR, DEV_EMAIL, NCBI_API_KEY

logger = logging.getLogger(__name__)

RETRIEVE_MAX = 10**7
EXTRA_QUERY_SLEEP_TIME = 2
RECORD_FIELDS = ['PMID', 'TI', 'AB', 'MH']


def validate_settings(email=DEV_EMAIL, api_key=NCBI_API_KEY):
    if email is None:
        logger.warning(f"Email is None! Is your .env file set correctly?")
        logger.info(f"Email adress used will be {Entrez.email} ")
        logger.warning(
            "Not setting a valid email adress can cause your IP to be blocked without ncbi having a way to contact you"
        )
    elif email == "your.devs.email@email.domain":
        logger.info(f"Email is still default {email}")
        logger.info(f"Email adress used will be {Entrez.email} ")
        logger.warning(
            "Not setting a valid email adress can cause your IP to be blocked without ncbi having a way to contact you"
        )
    else:
        Entrez.email = email
        logger.info(f"Email adress is set to {email}")

    logger.info(f"Max tries: {Entrez.max_tries}")
    logger.info(f"Sleep between tries: {Entrez.sleep_between_tries} s")

    if NCBI_API_KEY is None:
        logger.warning(f"NCBI_API_KEY is None! Is .env file set correctly?")
        logger.info("NO API KEY SET")
    elif NCBI_API_KEY != "s0meapikey":
        covered_key = api_key[:5] + ((len(api_key) - 10) * '*') + api_key[-5:]
        Entrez.api_key = api_key
        logger.info(f"NCBI API KEY set to: {covered_key}")
    else:
        logger.info("NO API KEY SET")


def _buildquery(keywords_subset,
                mesh_include=None,
                mesh_exclude=None,
                fromlist=True) -> str:

    keywords = ' OR '.join([
        f'({keyword.lower()})' for keyword in keywords_subset
    ]).replace('\n', ' ')

    include = ' AND '.join([f'{term.lower()}[mh]' for term in mesh_include
                            ]).replace('\n', ' ') if mesh_include else ""
    exclude = ' OR '.join([f'{term.lower()}[mh]' for term in mesh_exclude
                           ]).replace('\n', ' ') if mesh_exclude else ""

    query = f"({keywords})"
    query = f"{query} AND {include}" if include else query
    query = f"{query} NOT {exclude}" if exclude else query

    return query


def find_ids(
    search_keywords,
    # pre_formatted=False,
    mesh_include=[],
    mesh_exclude=[],
    fromlist=True,
    mindate='1500/01/01',  # the beginning of TIME
    email=DEV_EMAIL,
    api_key=NCBI_API_KEY,
    save_ids=True,  # TODO unimplemented
    save_ids_path=Path(DEFAULT_RECORDS_DIR / 'ids.txt').resolve(),
    save_queries=True,  # TODO unimplemented
    save_queries_path=Path(DEFAULT_RECORDS_DIR / 'queries_used.txt').resolve(),
    validate=True,
    new_only=False,  # TODO unimplemented
    chunk_size=100) -> Set[Bio.Entrez.Parser.StringElement]:
    """Search for article IDS on PubMed through Entrez.esearch
    

    Args:
        search_keywords (List[str]): list of keywords to look for
        mesh_include (str OR List[str], optional): list of MeSH terms to include. Defaults to [].
        mesh_exclude (str OR List[str], optional): list of MeSH terms to exclude. Defaults to [].
        fromlist (bool, optional): @notimplemented whether to build formatted string from MeSH terms or concatenate. Defaults to True.
        mindate (string, optional) The minimum date from which to search for articles. Will be overwritten with the current data if new_only==True. Defaults to '1500/01/01'
        email ([type], optional): email to be set into Bio.Entrez.email. Used to contact user when errors occur. Defaults to DEV_EMAIL.
        api_key ([type], optional): NCBI API key, get one at https://www.ncbi.nlm.nih.gov/account/ on the settings page . Defaults to NCBI_API_KEY, loaded from `.env` in taxmap.__init__ which will result in `validate_settings()` to not set an API key
        save_ids (bool, optional):  wheter to save the ids to a txtfile. Defaults to True.
        save_ids_path (str, optional):  the path to save ids to. Defaults to "./data/records/ids.txt"
        save_queries (bool, optional):  wheter to save the queries to a txtfile. Defaults to True.
        save_queries_path (str, optional): the path to save queries to. Defaults to "./data/records/queries_used.txt"
        new_only (bool, optional): . Defaults to False.
        chunk_size (int, optional) size of amount of keywords to use per query. Lower if each keyword is long on its own. Defaults to 100.
    Returns:
        the set of unique ids (of type Bio.Entrez.Parser.StringElement) found for all queries sent to pubmed

    Todo:
        Args:
            pre_formatted (bool, optional): Whether the list of keywords is preformatted, uses _buildquery by default to format each query chunk with MeSH terms. Defaults to False.


    """
    if validate:
        validate_settings(email, api_key)

    if fromlist:
        keywords = sorted(set(search_keywords))

    ids = set()
    queries_sent = list()
    for index in range(0, len(keywords), chunk_size):

        print(
            f'Looking for ids for queries {index: 6} to {index+chunk_size: 6} . . ',
            end=' ')

        chunk = keywords[index:index + chunk_size]

        query = _buildquery(chunk, mesh_include, mesh_exclude)

        handle = Entrez.esearch(db='pubmed',
                                term=query,
                                email=Entrez.email,
                                mindate=mindate,
                                retmax=RETRIEVE_MAX,
                                usehistory='y',
                                api_key=Entrez.api_key)
        results = Entrez.read(handle)
        handle.close()
        ids_found = results['IdList']
        ids = ids | set(ids_found)
        queries_sent.append(query)
        print(f'Found {len(ids_found): 8} ids')

    if save_ids:
        save_list(ids, save_ids_path)
        print(f'stored {len(ids)} ids to {save_ids_path}')
    if save_queries:
        save_list(queries_sent, save_queries_path)
        print(f'stored queries to {save_queries_path}')

    save_list(search_keywords,
              Path(DEFAULT_RECORDS_DIR / 'keywords_used.txt').resolve())

    return ids


def download_abstracts(
    ids,
    batch_size=5000,
    filename=Path(DEFAULT_RECORDS_DIR / 'records.txt').resolve(),
    errors_filename=Path(DEFAULT_RECORDS_DIR / 'errors.txt').resolve()):
    """Download all records from PUBMED using their corresponding ID

    Args:
        ids (List[int] or Str[int]): iterable type of Bio.Entrez.Parser.StringElement (!?)
        batch_size (int, optional): [description]. Defaults to 5000.
        filename ([type], optional): [description]. Defaults to Path(DEFAULT_RECORDS_DIR/'records.gzip').resolve().
    
    Todo: lots!
    """
    validate_settings()
    start_time = time.perf_counter()

    stats = Counter()
    batch_size = batch_size

    total_ids = len(ids)
    ids = sorted(ids)

    with open(filename, "w+") as output_file, open(errors_filename,
                                                   "w+") as errorfh:
        output_file.write("\t".join(RECORD_FIELDS) + "\n")

        for start in range(0, total_ids + 1, batch_size):

            batch = ids[start:start + batch_size]
            query = ",".join([str(Id) for Id in batch])

            batch_time_start = time.perf_counter()

            print(
                f"Downloading records {start+1: 8} to {start+batch_size if start+batch_size < total_ids else total_ids: 8} out of {total_ids}"
            )

            tries = 0
            fetch_handle = None
            while tries < 5:
                tries += 1
                try:
                    fetch_handle = Entrez.efetch(
                        db='pubmed',
                        id=query,
                        rettype="medline",  #todo check with 'abstract'
                        retmode="text",
                        retmax=batch_size * 2,
                        api_key=Entrez.api_key,
                        email=Entrez.email)
                    if fetch_handle:
                        break
                except HTTPError as e:
                    print(f'HTTPError {e.code}, Waiting a little extra...')
                    time.sleep(10)
                    
            records = Medline.parse(fetch_handle)

            for record in records:
                if not 'AB' in record.keys():
                    errorfh.write("\t".join([record['PMID'], "No abstract"]) +
                                  "\n")
                    stats["NO_ABSTRACT"] += 1
                    continue
                elif not 'TI' in record.keys():
                    record['TI'] = "No Title"
                    errorfh.write("\t".join([record['PMID'], "No Title"]) +
                                  "\n")
                    stats["NO_TITLE"] += 1

                if 'MH' in record.keys():
                    record['MH'] = ";;".join(record['MH'])
                else:
                    record['MH'] = 'NOMESHTERMS'
                output_file.write(
                    "\t".join([record[field]
                               for field in RECORD_FIELDS]) + '\n')
            fetch_handle.close()

            batch_time_end = time.perf_counter()
            duration = time.strftime(
                '%M minutes %S seconds',
                time.gmtime(batch_time_start - batch_time_end))
            print(
                f'. . . took {duration}. {total_ids-start: 8} records left to go'
            )

            time.sleep(EXTRA_QUERY_SLEEP_TIME)

    end_time = time.perf_counter()
    duration = time.strftime('%Dd %Hh %Mm %Sa ',
                             time.gmtime(start_time - end_time))
    print(f'It took {duration} to download {total_ids} records')
    for item, count in stats.most_common():
        print(item, count)
    print(f'The records can be found in {filename}')
