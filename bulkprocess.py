
import collections
import gzip
import logging
import os
import pickle
import sys
import time
import traceback
from collections import Counter, defaultdict
from urllib.error import HTTPError

import dotenv
import nltk
import spacy
import textacy
from Bio import Entrez, Medline

import bulkanalyser as ba
import bulkfilters as bf
from braintaxmap.config import dev_email
from braintaxmap.tools import load_verbs, timethisfunc
from bulkconstants import (BRAIN_FUNCTIONS, BRAIN_STRUCTURES, DATA_DIR,
                           MEDLINE_RESULTS_FILE)


def _medlineResultFhandle():
    with gzip.open(MEDLINE_RESULTS_FILE, "rb") as gzipfh:
        for line in gzipfh:
            yield str(line, 'utf-8')


def stored_records():
    return Medline.parse(_medlineResultFhandle())


def _read_amt_stored():
    records = stored_records()
    c = 0
    for r in records:
        c += 1

    return c


class MedlineAnalyser(object):
    threshold = 1  # min amt of occurences for stats

    # python -m spacy download en_core_web_lg 750Mb
    nlp = spacy.load("en_core_web_lg")
    verbs = load_verbs('data' + os.sep + '1000-verbs-set.txt')

    # 'Brain', 'Amygdala', 'Depression'
    mesh_include_list = ['Rats', 'Rodent', 'Mice', 'Muridae']
    mesh_exclude_list = ['Humans']
    relevant_terms = [set(BRAIN_FUNCTIONS), set(BRAIN_STRUCTURES)]

    def __init__(self):
        self._startup()
        self.finished = False

    def _set_initial_stats(self):
        self.global_stats = Counter()
        self.excluded_by = defaultdict(int)
        self.histogram_dicts = defaultdict(Counter)
        self.initial_stats = {
            "records_stored": _read_amt_stored(),
            # "Medline_file_size":os.path.getsize(MEDLINE_RESULTS_FILE)
        }

    def _print_initial_stats(self):
        print('INITIAL STATS:')
        for k, v in self.initial_stats.items():
            print(f'{k}:{v}')
        print(28*'-')

    def _startup(self):
        self._set_initial_stats()
        self._print_initial_stats()

    def _finish(self, reason='finished'):
        print(f'rounding up because of: {reason}')
        self.output_stats()
        self.finished = True
        print('done')

    @timethisfunc
    def run(self):
        try:
            self.run_analysis()
        except KeyboardInterrupt:
            self._finish('User stopped execution')

        if not self.finished:  # just in case, but hopefully unreachable
            print('finish missed. finishing up')
            self._finish()

    def run_analysis(self):
        print('starting to run analysis')
        records = stored_records()

        try:
            while True:
                record = next(records)
                self.analyse_record(record)
                self.global_stats['records analysed'] += 1
        except StopIteration:
            print('analysis done')
            self._finish()
            return

    @property
    def percentage_done(self):
        p = (self.global_stats['records analysed'] /
             self.initial_stats['records_stored']) * 100
        return p

    def analyse_record(self, record):
        """All things to be done with a record"""
        print(
            f"PMID:{record['PMID']:>10} ({self.percentage_done:>5.2f}%) filtering..", end='')
        if self._ApplyFiltersOnRecord(record):
            print()
            return

        print('.analysing abstract')
        self._analyse_abstract(record)

    def _ApplyFiltersOnRecord(self, record):
        """Returns TRUE if record does not pass any out of all requirements
        filters are like this
        'filtername' : (condition, on/off)
        """
        record_filters = [
            bf.NoAbstract(record, 1),
            bf.HasRelevantMeshterms(
                record, 0, self.mesh_include_list, self.mesh_exclude_list),
            bf.RelevantAbstract(
                record, 0, self.relevant_terms, force_one_of_each=False)
        ]

        if any([f.value for f in record_filters if f.value]):
            reasons = [f.__name__() for f in record_filters if f.value]
            # do something wiht pmid and reasons
            self._count_exclude_by(reasons)
            return True
        else:
            return False

    def _count_exclude_by(self, reasons, pmid=None):
        """TODO: see if +=pmid is viable (memorywise)"""
        for reason in reasons:
            self.excluded_by[reason] += 1

    def _analyse_abstract(self, record):

        for sentence in nltk.sent_tokenize(record['AB']):
            doc = self.nlp(sentence)

            for nounphrase in doc.noun_chunks:
                self.histogram_dicts['NOUNPHRASES'][nounphrase.text.lower(
                )] += 1

            for token in doc:
                ba.is_verb_fromListSpacy(
                    token, self.verbs, self.histogram_dicts)
                ba.is_verb_fromSpacy(token, self.histogram_dicts)

            for word in nltk.word_tokenize(sentence):
                ba.is_verb_fromListNLTK(word, self.verbs, self.histogram_dicts)

            ba.find_SVOs(sentence, self.nlp, self.histogram_dicts)

    def output_stats(self):

        for name, counter in self.histogram_dicts.items():
            with open(DATA_DIR+f'stats_{name}.txt', 'w+') as fh:
                fh.write(f"# STATS FOR {name.upper()} !\n\n")
                fh.write(
                    f"records analysed for this data: {self.global_stats['records analysed']} \n")
                for item, count in counter.most_common():
                    fh.write(f'{item}, {count}\n')


if __name__ == '__main__':
    m = MedlineAnalyser()
    m.run()
