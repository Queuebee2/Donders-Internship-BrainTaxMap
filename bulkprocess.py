
import collections
import gzip
import logging
import os
import pickle
import re
import shutil
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
from braintaxmap.tools import (dsm5parse, load_verbs, read_excluded,
                               read_included, readicd11simpleTabulation,
                               timethisfunc)
from bulkconstants import (BRAIN_FUNCTIONS, BRAIN_STRUCTURES, DATA_DIR,
                           MEDLINE_RESULTS_FILE, STATS_OUT_DIR)


def _medlineResultFhandle():
    with gzip.open(MEDLINE_RESULTS_FILE, "rb") as gzipfh:
        for line in gzipfh:
            yield str(line, 'utf-8')


def stored_records():
    return Medline.parse(_medlineResultFhandle())


def read_pmids_stored():
    records = stored_records()
    pmids = set()
    for r in records:
        pmids.add(r['PMID'])
    return pmids


def _read_amt_stored():
    records = stored_records()
    c = 0
    for r in records:
        c += 1

    return c


class MedlineAnalyser(object):
    threshold = 1  # min amt of occurences for stats (unsused)
    verbosity = 4  # 4 is max, 0 should be off

    def _setup_lists(self):

        # python -m spacy download en_core_web_lg 750Mb
        self.nlp = spacy.load("en_core_web_lg")
        # custom set of verbs
        self.verbs = load_verbs(os.path.join(DATA_DIR, '1000-verbs-set.txt')) - read_excluded() | read_included()
        # brainstructures of mouse http://help.brain-map.org/display/api/Downloading+an+Ontology%27s+Structure+Graph
        self.BRAIN_STRUCTURES=set(BRAIN_STRUCTURES)

        # https://brainmap.org/taxonomy/behaviors.html
        # (to be expanded, also look for 'disorders')
        self.BRAIN_FUNCTIONS=set(BRAIN_FUNCTIONS)

        # https://www.psychiatry.org/File%20Library/Psychiatrists/Practice/DSM/APA_DSM-5-Contents.pdf
        self.DISORDERS_DSM5=dsm5parse()
        # icd11 https://icd.who.int/en
        self.DISORDERS_ICD11=readicd11simpleTabulation()



        mesh_include_list=['Rats', 'Rodent', 'Mice', 'Muridae']
        mesh_exclude_list=['Humans']
        # 'Brain', 'Amygdala', 'Depression' etc
        relevant_terms=[set(BRAIN_FUNCTIONS), set(BRAIN_STRUCTURES)]

    def __init__(self):
        self._startup()
        self.finished=False

    def _set_initial_stats(self):
        self.global_stats=Counter()
        self.excluded_by=defaultdict(int)
        self.histogram_dicts=defaultdict(Counter)
        self.initial_stats={
            "records_stored": _read_amt_stored(),
            # "Medline_file_size":os.path.getsize(MEDLINE_RESULTS_FILE)
        }

    def _print_initial_stats(self):
        print('INITIAL STATS:')
        for k, v in self.initial_stats.items():
            print(f'{k}:{v}')
        print(28*'-')


    def _set_main_lists(self):
        self._setup_lists()
        self.HIT_TUPLES=list()
        self.DISORDERS=set()
        self.DISORDERS.update(self.DISORDERS_DSM5)
        for main_d, sub_d in self.DISORDERS_ICD11.items():
            self.DISORDERS.update(sub_d)
            self.DISORDERS.update({main_d})
        print(f'Merged disorders into one set. Length: {len(self.DISORDERS)}')

    def _create_regex(self, name, list):
        items=sorted(list, key=len)
        pattern_string=r"\b|".join([str(x) for x in items]) + r"\b"
        re_pattern=re.compile(pattern_string, re.IGNORECASE)
        self.REGEX_DICT[name]=re_pattern

        if self.verbosity > 2: print(
            f'compiled regex of length:{len(pattern_string)}, size={sys.getsizeof(re_pattern)} for {name}')

    def _set_regexes(self):
        self.REGEX_DICT=defaultdict(re.Pattern)

        self._create_regex('structures', self.BRAIN_STRUCTURES)
        self._create_regex('verbs', self.verbs)
        self._create_regex('functions', self.BRAIN_FUNCTIONS)
        self._create_regex('disorders', self.DISORDERS)


    def _startup(self):
        self._set_initial_stats()
        self._set_main_lists()
        self._set_regexes()
        self._refresh_output()
        self._print_initial_stats()

    def _finish(self, reason='finished'):
        print(f'rounding up because of: {reason}')
        self.output_stats()
        self.finished=True
        print('done')

    @ timethisfunc
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
        records=stored_records()

        try:
            while True:
                record=next(records)
                self.analyse_record(record)
                self.global_stats['records analysed'] += 1
        except StopIteration:
            print('analysis done')
            self._finish()
            return

    @ property
    def percentage_done(self):
        p=(self.global_stats['records analysed'] /
             self.initial_stats['records_stored']) * 100
        return p

    def analyse_record(self, record):
        """All things to be done with a record"""
        if self.verbosity > 3: print(
            f"PMID:{record['PMID']:>10} ({self.percentage_done:>5.2f}%) filtering..", end='')
        if self._ApplyFiltersOnRecord(record):
            return

        if self.verbosity > 3: print('.analysing abstract')
        self._analyse_abstract(record)

    def _ApplyFiltersOnRecord(self, record):
        """Returns TRUE if record does not pass any out of all requirements
        filters are like this
        'filtername' : (condition, on/off)
        """
        record_filters=[
            bf.NoAbstract(record, 1),
            bf.HasRelevantMeshterms(
                record, 0, self.mesh_include_list, self.mesh_exclude_list),

            # relevant abstract currently disabled in bulkfilters.py -> always returns True
            bf.RelevantAbstract(
                record, 0, self.relevant_terms, force_one_of_each=False)
        ]

        if any([f.value for f in record_filters if f.value]):
            reasons=[f.__name__() for f in record_filters if f.value]
            # do something wiht pmid and reasons
            self._count_exclude_by(reasons)
            return True
        else:
            return False

    def _count_exclude_by(self, reasons, pmid=None):
        """TODO: see if +=pmid is viable (memorywise)"""
        for reason in reasons:
            self.excluded_by[reason] += 1

    def __previous_methods(self, sentence):

            doc=self.nlp(sentence)
            ba.contains_icd11(sentence, doc, self.histogram_dicts)
            ba.contains_dsm5(sentence, doc, self.histogram_dicts)
            for nounphrase in doc.noun_chunks:
                self.histogram_dicts['NOUNPHRASES'][nounphrase.text.lower(
                )] += 1

            for token in doc:
                ba.is_verb_fromListSpacy(
                    token, self.verbs, self.histogram_dicts)
                ba.is_verb_fromSpacy(token, self.histogram_dicts)

            for word in nltk.word_tokenize(sentence):
                ba.is_verb_fromListNLTK(word, self.verbs, self.histogram_dicts)

            ba.find_SVOs(doc, self.histogram_dicts)

    def _current_methods(self, pmid, sentence):

        doc=self.nlp(sentence)

        ba.new_method(
            pmid,
            sentence,
            doc,
            self.REGEX_DICT,
            self.histogram_dicts,
            self.HIT_TUPLES)


    def _analyse_abstract(self, record):
        """Any functions here come from bulkanalyser.py
        toggle them off manually there
        """
        for sentence in nltk.sent_tokenize(record['AB']):

            self._current_methods(record['PMID'], sentence)

    def _refresh_output(self):
        try:
            shutil.copyfile(os.path.join(STATS_OUT_DIR, f'pmid hits.txt'), os.path.join(
                STATS_OUT_DIR, f'pmid hits - BACKUP.txt'))
        except:
            print('no backup made as no previous stats were stored')
        with open(os.path.join(STATS_OUT_DIR, f'pmid hits.txt'), 'w+'): pass
        print('refreshed outputfile:', os.path.join(
            STATS_OUT_DIR, f'pmid hits - BACKUP.txt'))

    def output_one_pmid_hit(self, hit):
        with open(os.path.join(STATS_OUT_DIR, f'pmid hits.txt'), 'a+') as fh:
            fh.write(f'{hit}\n')

    def _output_hitstrings(self):
        with open(os.path.join(STATS_OUT_DIR, f'pmid hits.txt'), 'a+') as fh:
            for pmid, struct, verb, disfunc_name in self.HIT_TUPLES:
                fh.write(f'{pmid};;; {struct};;; {verb};;; {disfunc_name}\n')

    def output_stats(self):
        print('outputting stats . . .')

        self._output_histograms()
        self._output_hitstrings()
        self._output_exclude_reasons()

        print('done outputting stats')

    def _output_exclude_reasons(self):
        stats_title_string=f"## STATS EXCLUDE REASONS!"
        headers_string='reason, count'
        total_records_string=f"## records analysed for this data: {self.global_stats['records analysed']}"
        total_excluded=sum(val for val in self.excluded_by.values())
        total_excluded_string=f"## records excluded for this data: {total_excluded}"
        print(stats_title_string)
        print(total_records_string)
        with open(os.path.join(STATS_OUT_DIR, f'Exclude reasons counts'), 'w+') as fh:
            fh.write(f"{stats_title_string}\n\n")
            fh.write(f"{total_records_string}\n")
            fh.write(f"{total_excluded_string}\n\n")
            fh.write(f"{headers_string}\n")
            for reason, count in self.excluded_by.items():
                fh.write(f'{reason}, {count}\n')

    def _output_histograms(self):
        for name, counter in self.histogram_dicts.items():
            stats_title_string=f"## STATS FOR {name.upper()} !"
            headers_string='item, occurences'
            total_records_string=f"## records analysed for this data: {self.global_stats['records analysed']}"
            print(stats_title_string)
            print(total_records_string)

            with open(os.path.join(STATS_OUT_DIR, f'stats_{name}.txt'), 'w+') as fh:
                fh.write(f"{stats_title_string}\n\n")
                fh.write(f"{total_records_string}\n\n")
                fh.write(f"{headers_string}\n")
                for item, count in counter.most_common():
                    fh.write(f'{item}, {count}\n')

if __name__ == '__main__':
    m=MedlineAnalyser()
    m.run()
