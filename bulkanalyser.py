
from collections import defaultdict
import nltk
import spacy
import re
import textacy
from braintaxmap.tools import readcd11simpleTabulation, dsm5parse
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


cd11_dict = readcd11simpleTabulation()
def contains_cd11(sentence, sentence_doc, stats_dict):
    sentence=sentence.lower()
    for k, s in cd11_dict.items():
        if k.lower() in sentence:
            svo_triples = textacy.extract.subject_verb_object_triples(sentence_doc)
            svo_triples = list(svo_triples)
            stats_dict['foundcd11'][f'{k};;; {sentence};;; SVOs: {svo_triples}'] += 1
        for v in s:
            if v.lower() in sentence:
                svo_triples = textacy.extract.subject_verb_object_triples(sentence_doc)
                svo_triples = list(svo_triples)
                stats_dict['foundcd11'][f'{v};;; {sentence};;; SVOs: {svo_triples}'] += 1

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

def new_method(pmid, sentence, doc, structures, verbs, functions, disorders, stats_dict):
    hitstrings= []
    hits=defaultdict(list)

    # run all regexes against the sentence. This takes a ridiculous amount of time.
    for name, l in zip(['structures','verbs','functions','disorders'],
                          [structures, verbs, functions, disorders]):
        for word in l:
            match = re.search(rf"{word}\b", sentence, re.IGNORECASE)
            if match:
                hits[name].append(match.group(0))

    # check if at least a structure and verb an one of disorder/function is present
    if len(hits['structures']) > 0 and  len(hits['verbs']) > 0 and (
        len(hits['functions']) > 0 or len(hits['disorders']) > 0):

        # this creates all combinations of found items for this sentence
        # its because the SVO detector does not work. Hits should be manually reviewed
        for s in hits['structures']:
            for v in hits['verbs']:
                for a, b in zip(hits['disorders'], hits['functions']):
                    # hitstrings will be returned and immediately stored to a file
                    # todo: just keep them in class? because they are very few.
                    if a:
                        hitstrings.append(f'{pmid};;; {s};;; {v};;; {a}')
                    if b:
                        hitstrings.append(f'{pmid};;; {s};;; {v};;; {b}')
        # keep track of items we found
        for k, v in hits.items():
            for hit in v:
                stats_dict[f'{k} found'][hit]+=1
    else:   
        # keep track of items we found, but arent used
        for k, v in hits.items():
            for hit in v:
                stats_dict[f'{k} unused'][hit]+=1
    
    return hitstrings