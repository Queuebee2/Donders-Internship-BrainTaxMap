
import gzip
import os
import sys
import time
from pickle import load
from urllib.error import HTTPError

from Bio import Entrez
from dotenv import find_dotenv, load_dotenv

from braintaxmap.tools import timethisfunc_dhms
from braintaxmap.config import dev_email
# structures and behaviour funcs
from bulkconstants import (API_KEY, BRAIN_FUNCTIONS, BRAIN_STRUCTURES,
                           DATA_DIR, MEDLINE_RESULTS_FILE, STATS_OUT_DIR)
from bulkprocess import read_pmids_stored

print(
    f'loaded {len(BRAIN_STRUCTURES)} structures and {len(BRAIN_FUNCTIONS)} behaviours')

# script has been used, using this to prevent accidentally overwriting stuff.
print(' OR '.join(
    [f'({query})' for query in BRAIN_STRUCTURES+BRAIN_FUNCTIONS]).replace('\n', ' '))


# Entrez settings

# private API KEY

print('initial email:', Entrez.email, 'initial max tries:', Entrez.max_tries)
Entrez.email = dev_email
Entrez.max_tries = 15
Entrez.sleep_between_tries = 20
Entrez.api_key = API_KEY
print('email set to', Entrez.email, 'max tries set to', Entrez.max_tries)

retmax = 10**7
EXTRA_QUERY_SLEEP_TIME = 2  # seconds
CONTINUE_FROM_PREV = True

QUERIES = BRAIN_STRUCTURES

with open(os.path.join(*[DATA_DIR, 'pruned_structures_list.txt']), 'r') as fh:
    for line in fh:
        QUERIES.append(line.strip())


def _save_last_pmid(pmid):
    with open('last-pmid.txt', 'w+') as fh:
        fh.write(f'{pmid}')


def _get_last_pmid():
    with open('last-pmid.txt', 'r+') as fh:
        pmid = fh.read()
        return pmid


def _load_list_from_file(filename):
    result = list()
    with open(filename, 'r+') as fh:
        for line in fh:
            result.append(line.strip())
    return result


def _store_list_in_file(lst, filename):
    with open(filename, 'w+') as fh:
        for item in lst:
            fh.write(f'{item}\n')


def _get_stored_pmids():
    try:
        return _load_list_from_file(os.path.join(*[DATA_DIR, 'unique-pmids-found.txt']))
    except FileNotFoundError:
        return []


def _store_pmids(pmids):
    _store_list_in_file(pmids, os.path.join(
        *[DATA_DIR, 'unique-pmids-found.txt']))


def get_all_ids(queries_list):

    ids = _get_stored_pmids()
    if len(ids) > 0:
        print(f'successfully loaded {len(ids)} ids from storage')
        return ids

    queries = sorted(list(set(queries_list)))

    chunksize = 100
    ids = set()
    for qstart in range(0, len(queries), chunksize):

        print(
            f'Looking for ids for queries {qstart: 6} to {qstart+chunksize: 6} . . ', end=' ')
        chunk = queries[qstart:qstart+chunksize]
        query = ' OR '.join(
            [f'({query})' for query in chunk]).replace('\n', ' ')

        handle = Entrez.esearch(
            db='pubmed',
            term=query,
            email=dev_email,
            mindate='1500/01/01',
            retmax=retmax,
            usehistory='y',
            api_key=API_KEY)

        result = Entrez.read(handle)
        handle.close()
        ids_found = result['IdList']
        ids = ids | set(ids_found)
        print(f'Found {len(ids_found): 8} ids')
        time.sleep(EXTRA_QUERY_SLEEP_TIME)

    _store_pmids(ids)
    print(f'stored {len(ids)} ids to storage')

    return list(ids)


# timing
tstart = time.perf_counter()

# get ids
print(f"Looking for {len(QUERIES)} queries might take approx. {EXTRA_QUERY_SLEEP_TIME*len(QUERIES) + len(QUERIES)} seconds")
ids = get_all_ids(QUERIES)
print(f"Found {len(ids)} ids")

prev_ids = read_pmids_stored()
print(
    f"Ids currently already stored: {len(prev_ids)} of which {len(prev_ids.intersection(ids))} overlap")

# TODO, previously stored does not matter unless the records are
# stored apart. All other IDs should be removed from the file or stored somewhere else ?

ids = set(ids) - prev_ids
print(f"There are {len(ids)} ids leftover after subtracting the stored ids")

print('max size =', sys.maxsize)
print('size of ids-set =', sys.getsizeof(ids))
print(f'found {len(ids)} ids for {len(QUERIES)} different queries')

# timing
tend = time.perf_counter()
duration = time.strftime('%Dd %Hh %Mm %Sa ', time.gmtime(tend-tstart))
print(f'That took {duration}')
tstart = time.perf_counter()  # timing next

# download medline records as text
print(f'Starting to download {len(ids)} items...')

failed_records = []
batch_size = 5000
print(f'Batch size is {batch_size}')
TOTAL_IDS = len(ids)
ids = sorted(list(ids))
with gzip.open(MEDLINE_RESULTS_FILE, "ab") as gzipfh:

    for start in range(0, len(ids)+1, batch_size):

        ids_chunk = ids[start:start+batch_size]
        query = ",".join([str(i) for i in ids_chunk])

        chunk_time_start = time.perf_counter()

        print(
            f"Downloading records {start+1} to {start+batch_size} out of {TOTAL_IDS}")

        tries = 0
        fetch_handle = None
        while tries < 5:
            tries += 1
            try:
                fetch_handle = Entrez.efetch(
                    db="pubmed",
                    id=query,
                    rettype="medline",
                    retmode="text",
                    retmax=batch_size*2,
                    api_key=API_KEY,
                    email=dev_email
                )
                if fetch_handle:
                    break
            except HTTPError as e:
                print(f'HTTPError {e.code}, Waiting a little extra...')
                time.sleep(15)

        # retrieve and write data
        data = fetch_handle.read()
        fetch_handle.close()
        gzipfh.write(data.encode())

        # time and stats
        chunk_time_end = time.perf_counter()
        duration = time.strftime('%M minutes %S seconds', time.gmtime(
            chunk_time_end-chunk_time_start))
        print(f'. . . took {duration}. {TOTAL_IDS-start: 8} ids left to go!')

        time.sleep(EXTRA_QUERY_SLEEP_TIME)


# timing
tend = time.perf_counter()
duration = time.strftime('%Dd %Hh %Mm %Sa ', time.gmtime(tend-tstart))
print(f'That took {duration} seconds for {TOTAL_IDS} records')

print('done')
