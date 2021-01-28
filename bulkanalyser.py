
import re
from collections import defaultdict
import sys
import nltk
import spacy
import textacy

from braintaxmap.tools import dsm5parse, readicd11simpleTabulation

"""Anything that isnt used can be disabled by just placing an empty return at the start"""


def find_SVOs(sentence_doc, stats_dict):
    return
    """Try to extract SVO from a sentence and return a list of found SVO triples (s,v,o)"""
    sentence_doc
    svo_triples = textacy.extract.subject_verb_object_triples(sentence_doc)
    svo_triples = list(svo_triples)
    for svo in svo_triples:
        keysvo = str(svo).lower()
        stats_dict['SVOs'][keysvo] += 1


def is_verb_fromListSpacy(word, verbs, stats_dict):
    return
    word = str(word).lower()
    if word in verbs:
        stats_dict['VERBS (spacy tokenized)'][word] += 1


def is_verb_fromListNLTK(word, verbs, stats_dict):
    return
    word = str(word).lower()
    if word.lower() in verbs:
        stats_dict['VERBS (nltk tokenized)'][word] += 1


def is_verb_fromSpacy(word, stats_dict):
    return
    if word.pos_ == 'VERB':
        stats_dict['VERBS (spacy token.pos_)'][word] += 1


def is_nounphrase(nounphrase, stats_dict):
    return
    stats_dict['nounphrases'][nounphrase.text.lower()] += 1


icd11_dict = readicd11simpleTabulation()
def contains_icd11(sentence, sentence_doc, stats_dict):
    sentence=sentence.lower()
    for k, s in icd11_dict.items():
        if k.lower() in sentence:
            svo_triples = textacy.extract.subject_verb_object_triples(sentence_doc)
            svo_triples = list(svo_triples)
            stats_dict['foundicd11'][f'{k};;; {sentence};;; SVOs: {svo_triples}'] += 1
        for v in s:
            if v.lower() in sentence:
                svo_triples = textacy.extract.subject_verb_object_triples(sentence_doc)
                svo_triples = list(svo_triples)
                stats_dict['foundicd11'][f'{v};;; {sentence};;; SVOs: {svo_triples}'] += 1

dsm5terms = dsm5parse()
def contains_dsm5(sentence, sentence_doc, stats_dict):
    sentence=sentence.lower()
    for t in dsm5terms:
        if t.lower() in sentence:
            svo_triples = textacy.extract.subject_verb_object_triples(sentence_doc)
            svo_triples = list(svo_triples)
            stats_dict['dsm5found'][f'{t};;; {sentence};;; SVOs: {svo_triples}'] += 1


def prev_method(sentence, doc, structures, verbs, functions, disorders, stats_dict):
    """FAILED because too slow"""

    hits=defaultdict(list)
    for name, l in zip(['structures','verbs','functions','disorders'],
                          [structures, verbs, functions, disorders]):
        for word in l:
            match = re.search(rf"{word}\b", sentence, re.IGNORECASE)
            if match:
                hits[name].append(match.group(0))

    if len(hits['structures']) > 0 and  len(hits['verbs']) > 0 and (
        len(hits['functions']) > 0 or len(hits['disorders']) > 0):
        svo_triples = textacy.extract.subject_verb_object_triples(doc)
        svo_triples = list(svo_triples)
        for svo in svo_triples:
            svo_as_string = str(svo)
            for s in hits['structures']:
                for v in hits['verbs']:
                    for i in zip(hits['disorders'], hits['functions']):
                        if (re.search(rf"{s}\b", svo_as_string, re.IGNORECASE) and 
                            re.search(rf"{v}\b", svo_as_string, re.IGNORECASE) and
                            re.search(rf"{i}\b", svo_as_string, re.IGNORECASE)):
                            stats_dict['SVOs'][svo_as_string.lower()] += 1
    
        # keep track of items we found
        for k, v in hits.items():
            for hit in v:
                stats_dict[f'{k} found'][hit]+=1

    else:   
        # keep track of items we found, but arent used
        for k, v in hits.items():
            for hit in v:
                stats_dict[f'{k} unused'][hit]+=1

def new_method(pmid, sentence, regex_dict, stats_dict, HIT_TUPLES, verbose=False):

    hits=defaultdict(list)
    for regname, regexpattern in regex_dict.items():
        re_hit = regexpattern.search(sentence)
        if re_hit:
            hits[regname].append(re_hit.group(0))
    
    # check if at least a structure and verb an one of disorder/function is present
    if len(hits['structures']) > 0 and  len(hits['verbs']) > 0 and (
        len(hits['functions']) > 0 or len(hits['disorders']) > 0):

        # this creates all combinations of found items for this sentence
        # its because the SVO detector does not work. Hits should be manually reviewed
        for struct in hits['structures']:
            for verb in hits['verbs']:
                # hitstrings will be returned and immediately stored to a file
                # todo: just keep them in class? because they are very few.
                for disorder_name in hits['disorders']:
                    HIT_TUPLES['disorders'].append((pmid,struct,verb,disorder_name))
                    if verbose: print('found tuple.. ', (pmid,struct,verb,disorder_name))
                for func_name in  hits['functions']:
                    HIT_TUPLES['functions'].append((pmid,struct,verb,func_name))
                    if verbose: print('found tuple.. ', (pmid,struct,verb,func_name))

        # keep track of items we found
        for k, verb in hits.items():
            for hit in verb:
                stats_dict[f'{k} found'][hit]+=1
    else:   
        # keep track of items we found, but arent used
        for k, verb in hits.items():
            for hit in verb:
                stats_dict[f'{k} unused'][hit]+=1
    
