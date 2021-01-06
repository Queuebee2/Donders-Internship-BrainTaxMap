
from pickle import load
from Bio import Entrez, Medline
from braintaxmap.config import dev_email
import time
import nltk
import sys
import os
import gzip
from bulkconstants import DATA_DIR, MEDLINE_RESULTS_FILE



# structures and behaviour funcs
from braintaxmap.data_prep import load_brainstructures, load_behaviours
structures_pickle_path = r'data\\extracted-json-nodes.pickle'
behaviours_path = r'data\\list_of_behaviour_hirarchy.txt'
brainstructures = [n.lower() for n in  load_brainstructures(structures_pickle_path)]
# brainbehaviours = [n.lower() for n in load_behaviours(behaviours_path)]
brainstructures.remove('root')

# Entrez settings

## private API KEY
API_KEY = open(r"..\\NCBI-API-KEY.txt").read()[:-1]


print('initial email:', Entrez.email, 'initial max tries:', Entrez.max_tries)
Entrez.email = dev_email
Entrez.max_tries = 15
Entrez.sleep_between_tries = 20
Entrez.api_key = API_KEY
print('email set to', Entrez.email, 'max tries set to', Entrez.max_tries)


retmax=10**7
QUERY_INTERVAL = 2 #seconds
CONTINUE_FROM_PREV=True

queries = [
    'barrel cortex', 
    '(amygdala) AND (depression)',
    '(barrel cortex) AND (depression)',
    '(uncertainty) and (brain) and (behaviour)',
    '(language) and (brain)' ,
    '(EEG) and (behaviour)',
    '(Brain) and (disease) and (behaviour)'
    ] + brainstructures 

def update_ids(all_ids):
    ids_in_storage = set()
    with open(MEDLINE_RESULTS_FILE, 'rt') as fh:
        for line in fh:
            # example line PMID- 20033235
            if line.startswith("PMID-"):
                pmid = int(line.strip().split(' ')[1])
                ids_in_storage.add(pmid)
    
    return all_ids - ids_in_storage

def get_all_ids(queries):
    
    chunksize=100
    ids = set()
    for qstart in range(0, len(queries), chunksize):
        chunk = queries[qstart:qstart+chunksize]
        query=' OR '.join([f'({query})' for query in chunk]).replace('\n', ' ')

        handle = Entrez.esearch(
                db='pubmed',
                term=query,
                email=dev_email,
                mindate='1500/01/01',
                retmax=retmax,
                usehistory='y',
                api_key=API_KEY)
            
        result = Entrez.read(handle)
        count = int(results['Count'])
        webenv = results["WebEnv"]
        query_key = results["QueryKey"]
    
    return result

# timing
tstart = time.perf_counter()

# get ids
print(f"Looking for {len(queries)} queries might take approx. {QUERY_INTERVAL*len(queries) + len(queries)*2} seconds")
results = get_all_ids(queries)
count = int(results['Count'])
webenv = results["WebEnv"]
query_key = results["QueryKey"]

ids = set(results['IdList'])
print(f'count == len(ids-list)? {len(ids) == count}')
print('max size =',sys.maxsize)
print('size of ids set =',sys.getsizeof(ids))
print(f'found {len(ids)} ids for {len(queries)} different queries')

# timing
tend = time.perf_counter()
duration = tend-tstart
print(f'That took {duration} s')
tstart = time.perf_counter() # timing next 

# download medline records as text
print('Starting downloading items...')

failed_records = []
batch_size = 500

with gzip.open(MEDLINE_RESULTS_FILE, "wb") as gzipfh:
    for start in range(0, 500, batch_size):
        end = min(count, start + batch_size)
        print(f"Downloading records {start+1} to {end} out of {count}")
        fetch_handle = Entrez.efetch(
            db="pubmed",
            rettype="medline",
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key
        )
        data = fetch_handle.read()
        gzipfh.write(data.encode())


# timing
tend = time.perf_counter()
print(f'That took {duration} seconds for {count} records')



print('done')