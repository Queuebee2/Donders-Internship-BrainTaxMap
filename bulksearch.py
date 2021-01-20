
from braintaxmap.data_prep import load_previous
import gzip
import os
import sys
import time
from pickle import load
from urllib.error import HTTPError

import nltk
from Bio import Entrez, Medline

from braintaxmap.config import dev_email
# structures and behaviour funcs
from bulkconstants import (BRAIN_FUNCTIONS, BRAIN_STRUCTURES, STATS_OUT_DIR,
                           MEDLINE_RESULTS_FILE)
from bulkprocess import read_pmids_stored

print(
    f'loaded {len(BRAIN_STRUCTURES)} structures and {len(BRAIN_FUNCTIONS)} behaviours')

# script has been used, using this to prevent accidentally overwriting stuff.
print(' OR '.join(
    [f'({query})' for query in BRAIN_STRUCTURES+BRAIN_FUNCTIONS]).replace('\n', ' '))


# Entrez settings

# private API KEY
API_KEY = open(r"..\\NCBI-API-KEY.txt").read()[:-1]

print('initial email:', Entrez.email, 'initial max tries:', Entrez.max_tries)
Entrez.email = dev_email
Entrez.max_tries = 15
Entrez.sleep_between_tries = 20
Entrez.api_key = API_KEY
print('email set to', Entrez.email, 'max tries set to', Entrez.max_tries)

retmax = 10**7
EXTRA_QUERY_SLEEP_TIME = 2  # seconds
CONTINUE_FROM_PREV = True

queries = [
    'barrel cortex',
    '(amygdala) AND (depression)',
    '(barrel cortex) AND (depression)',
    '(uncertainty) and (brain) and (behaviour)',
    '(language) and (brain)',
    '(EEG) and (behaviour)',
    '(Brain) and (disease) and (behaviour)'
] + BRAIN_STRUCTURES

def _save_last_pmid(pmid):
    with open('last-pmid.txt', 'w+') as fh:
        fh.write(f'{pmid}')

def _get_last_pmid():
    with open('last-pmid.txt', 'r+') as fh:
        pmid = fh.read()
        return pmid

def get_all_ids(queries):
    
    chunksize = 100
    ids = set()
    for qstart in range(0, len(queries), chunksize):

        print(
            f'Looking for ids for queries {qstart} to {qstart+chunksize} . . ', end=' ')
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
        print(f'Found {len(ids_found)} ids')
        time.sleep(EXTRA_QUERY_SLEEP_TIME)

    return list(ids)


# timing
tstart = time.perf_counter()

# get ids
print(f"Looking for {len(queries)} queries might take approx. {EXTRA_QUERY_SLEEP_TIME*len(queries) + len(queries)} seconds")
ids = get_all_ids(queries)

prev_ids = read_pmids_stored()
print(f"Ids currently already stored: {len(prev_ids)} of which {len(prev_ids.intersection(ids))} overlap")
# TODO, previously stored does not matter unless the records are
# stored apart. All other IDs should be removed from the file or stored somewhere else ?

ids = set(ids) - prev_ids
print('quitting..')
quit()

print('max size =', sys.maxsize)
print('size of ids-set =', sys.getsizeof(ids))
print(f'found {len(ids)} ids for {len(queries)} different queries')

# timing
tend = time.perf_counter()
duration = time.strftime('%M minutes %S seconds', time.gmtime(tend-tstart))
print(f'That took {duration}')
tstart = time.perf_counter()  # timing next

# download medline records as text
print(f'Starting to download {len(ids)} items...')

failed_records = []
batch_size = 5000
print(f'Batch size is {batch_size}')
TOTAL_IDS = len(ids)
with gzip.open(MEDLINE_RESULTS_FILE, "wb") as gzipfh:

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
        print(f'. . . took {duration}. {TOTAL_IDS-start} ids left to go!')

        time.sleep(EXTRA_QUERY_SLEEP_TIME)


# timing
tend = time.perf_counter()
print(f'That took {duration} seconds for {TOTAL_IDS} records')

print('done')
