

import requests

from Bio.Entrez import esearch
import xml.etree.ElementTree as ET
import time
from dotenv import find_dotenv
import dotenv

IDS_SAVEFILE = 'ids-found.txt'
XML_SAVEFILE = 'xml-found.xml'
E_UTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def _maketerm(query):
    return query.replace(' ','+')

def _create_uri(tool_link, query, db, retmax=10**7, usehistory='y'):
    term = _maketerm(query)
    uri = f"{E_UTILS}{tool_link}?db={db}&term={term}&usehistory={usehistory}&retmax={retmax}"
    return uri



def _check_idsfile(ids_savefile=IDS_SAVEFILE):
    with open(ids_savefile, 'r') as fh:
        c=0
        for line in fh:
            c+=1
        return c

def esearch(query, db='pubmed', ids_savefile = IDS_SAVEFILE):
    """
    default db = pubmed
    """
    
    tool_link = 'esearch.fcgi'
    uri = _create_uri(tool_link, query, db)
    
    print(f'Looking for {query} with {uri}')

    response = requests.get(uri, stream=True)
    # if the server sent a Gzip or Deflate compressed response, decompress
    # as we read the raw stream:
    response.raw.decode_content = True

    events = ET.iterparse(response.raw)
    ids_found=0


    with open(ids_savefile, 'w') as fh:
        for event, elem in events:
            # do something with `elem`
            if elem.tag == 'Id':
                fh.write(f"{elem.text}\n")
                ids_found+=1
            if 'webenv' in str(elem.tag).lower():
                web_env = elem.text
                print(elem.tag)
                print(elem.text)
                dotenv.set_key(find_dotenv(), "WEB_ENV", web_env)
            if 'querykey' in str(elem.tag).lower():
                print(elem.tag)
                print(elem.text)
                query_key = elem.text
                dotenv.set_key(find_dotenv(), "QUERY_KEY", query_key)

    

    print(f'ids found for {query}: {ids_found} and stored in {ids_savefile}')
    ids_stored = _check_idsfile(ids_savefile)
    print(f'lines (ids) verified as stored in idsfile: {ids_stored} == {ids_found} ? {ids_found==ids_stored}')
    assert ids_found==ids_stored, "unequal amount of ids found and stored"

def _clear_previous_xml(xml_file):
    with open(xml_file, 'w+') as f:
        pass

def efetch( chunksize=25, idsfile=IDS_SAVEFILE):

    tool_link = 'efetch.fcgi?'
    
    idsread=0
    idslimit=_check_idsfile(idsfile)
    web_env = dotenv.get_key(find_dotenv(), 'WEB_ENV')
    query_key = dotenv.get_key(find_dotenv(), 'QUERY_KEY')
    _clear_previous_xml(XML_SAVEFILE)
    while idsread< idslimit:

        current_ids=[]
        
        with open(idsfile, 'r') as fh:
            for i, line in enumerate(fh):
                if i >= idsread:
                    id = line.strip()
                    current_ids+=[id]
                    if len(current_ids)>= chunksize:
                        idsread+=len(current_ids)
                        break
        
        ids = ",".join(current_ids)
        
        uri = f"{E_UTILS}{tool_link}db=pubmed&ids={ids}&rettype=medline&retmode=text&WebEnv={web_env}&query_key={query_key}" #&usehistory=y"

        print(f'looking for {len(current_ids)} ids from {idsread} to {idsread+len(current_ids)}')
        response = requests.get(uri, stream=True)
        print('got response . . . ')

        # # if the server sent a Gzip or Deflate compressed response, decompress
        # # as we read the raw stream:
        # response.raw.decode_content = True

        # events = ET.iterparse(response.raw)
        # ids_found=0
        # tags_found=set()
        # with open(XML_SAVEFILE, 'a+') as fh:
        #     for event, elem in events:
        #         # do something with `elem`
        #         # if elem.tag == '':
        #             # fh.write(f"{elem.text}\n")
        #             # ids_found+=1
        #         tags_found.add(elem.tag)
        #         if elem.tag == 'Title':
        #             print(elem.text)
        #         if elem.tag == '' 

        for line in response.iter_lines(delimiter=b'\n'):
            print(line)

        print('tags found:', tags_found)
        print(uri)
        quit()

    pass

esearch('brain')
efetch()