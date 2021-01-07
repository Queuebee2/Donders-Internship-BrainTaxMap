
import nltk
import spacy
import textacy


def find_SVOs(sentence, nlp, stats_dict):
    """Try to extract SVO from a sentence and return a list of found SVO triples (s,v,o)"""
    doc = nlp(sentence)
    svo_triples = textacy.extract.subject_verb_object_triples(doc)
    svo_triples = list(svo_triples)
    for svo in svo_triples:
        keysvo = str(svo).lower()
        stats_dict['SVOs'][keysvo] += 1


def is_verb_fromListSpacy(word, verbs, stats_dict):
    word = str(word).lower()
    if word in verbs:
        stats_dict['VERBS (spacy tokenized)'][word] += 1


def is_verb_fromListNLTK(word, verbs, stats_dict):
    word = str(word).lower()
    if word.lower() in verbs:
        stats_dict['VERBS (nltk tokenized)'][word] += 1


def is_verb_fromSpacy(word, stats_dict):
    if word.pos_ == 'VERB':
        stats_dict['VERBS (spacy token.pos_)'][word] += 1


def is_nounphrase(nounphrase, stats_dict):
    stats_dict['nounphrases'][nounphrase.text.lower()] += 1
