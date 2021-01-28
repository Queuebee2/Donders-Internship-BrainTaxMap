
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
                               timethisfunc_hms, timethisfunc_dhms)
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
    print_interval = 5000

    def _setup_lists(self):

        # python -m spacy download en_core_web_lg 750Mb
        self.nlp = None # spacy.load("en_core_web_lg")
        print('NLP disabled' if self.nlp is None else 'NLP enabled')
        # custom set of verbs

        self.verbs_excluded = read_excluded()
        self.verbs_included = read_included()
        self.verbs = load_verbs(os.path.join(
            DATA_DIR, '1000-verbs-set.txt')) - self.verbs_excluded | self.verbs_included
        # brainstructures of mouse http://help.brain-map.org/display/api/Downloading+an+Ontology%27s+Structure+Graph
        self.BRAIN_STRUCTURES = set(BRAIN_STRUCTURES)

        # https://brainmap.org/taxonomy/behaviors.html
        # (to be expanded, also look for 'disorders')
        self.BRAIN_FUNCTIONS = set(BRAIN_FUNCTIONS)

        # https://www.psychiatry.org/File%20Library/Psychiatrists/Practice/DSM/APA_DSM-5-Contents.pdf
        self.DISORDERS_DSM5 = dsm5parse()
        # icd11 https://icd.who.int/en
        self.DISORDERS_ICD11 = readicd11simpleTabulation()

        self.mesh_include_list = ['Rats', 'Rodent', 'Mice', 'Muridae']
        self.mesh_exclude_list = ['Humans']
        # 'Brain', 'Amygdala', 'Depression' etc
        self.relevant_terms = [set(BRAIN_FUNCTIONS), set(BRAIN_STRUCTURES)]

    def __init__(self):
        self._startup()
        self.finished = False

    def _set_initial_stats(self):
        self.global_stats = Counter()
        self.excluded_by_counts = defaultdict(int)
        self.excluded_by_pmids = defaultdict(list)
        self.histogram_dicts = defaultdict(Counter)
        self.initial_stats = {
            "records_stored": _read_amt_stored(),
            "print interval": self.print_interval,
            "verbosity level": self.verbosity,      
            # "Medline_file_size":os.path.getsize(MEDLINE_RESULTS_FILE)
        }

    def _print_initial_stats(self):
        print('INITIAL STATS:')
        for k, v in self.initial_stats.items():
            print(f'{k}:{v}')
        print(28*'-')

    def _set_main_lists(self):
        self._setup_lists()
        self.HIT_TUPLES = defaultdict(list)
        self.HIT_TUPLES.update({k: list() for k in ['disorders', 'functions']})
        self.DISORDERS = set()
        self.DISORDERS.update(self.DISORDERS_DSM5)
        for main_d, sub_d in self.DISORDERS_ICD11.items():
            self.DISORDERS.update(sub_d)
            self.DISORDERS.update({main_d})
        print(f'Merged disorders into one set. Length: {len(self.DISORDERS)}')

    def _create_regex(self, name, list):
        items = sorted(list, key=len)
        pattern_string = r"\b|".join([str(x) for x in items]) + r"\b"
        re_pattern = re.compile(pattern_string, re.IGNORECASE)
        self.REGEX_DICT[name] = re_pattern

        if self.verbosity > 2:
            print(
                f'compiled regex of length:{len(pattern_string)}, size={sys.getsizeof(re_pattern)} for {name}')

    def _set_regexes(self):
        self.REGEX_DICT = defaultdict(re.Pattern)

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
        self.finished = True
        print('done')

    @ timethisfunc_dhms
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

    @ property
    def print_allowed_by_interval(self):
        return self.global_stats['records analysed'] % self.print_interval == 0

    @ property
    def percentage_done(self):
        p = (self.global_stats['records analysed'] /
             self.initial_stats['records_stored']) * 100
        return p

    def analyse_record(self, record):
        """All things to be done with a record"""
        if self.verbosity > 3 and self.print_allowed_by_interval:
            print(
                f"PMID:{record['PMID']:>10} ({self.percentage_done:>5.2f}%) filtering..", end='')
        if self._ApplyFiltersOnRecord(record):
            return

        if self.verbosity > 3 and self.print_allowed_by_interval:
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
                record, 1, self.mesh_include_list, self.mesh_exclude_list),

            # relevant abstract currently disabled in bulkfilters.py -> always returns True
            bf.RelevantAbstract(
                record, 0, self.relevant_terms, force_one_of_each=False)
        ]

        if any([f.value for f in record_filters if f.value]):
            reasons = [f.__name__() for f in record_filters if f.value]
            # do something wiht pmid and reasons
            self._count_exclude_by(reasons, record['PMID'])

            return True
        else:
            return False

    def _count_exclude_by(self, reasons, pmid=None):
        """TODO: see if +=pmid is viable (memorywise)"""
        for reason in reasons:
            self.excluded_by_counts[reason] += 1
            self.excluded_by_pmids[reason] += [pmid]

    def __previous_methods(self, sentence):

        # doc = self.nlp(sentence)
        # ba.contains_icd11(sentence, doc, self.histogram_dicts)
        # ba.contains_dsm5(sentence, doc, self.histogram_dicts)
        # for nounphrase in doc.noun_chunks:
        #     self.histogram_dicts['NOUNPHRASES'][nounphrase.text.lower(
        #     )] += 1

        # for token in doc:
        #     ba.is_verb_fromListSpacy(
        #         token, self.verbs, self.histogram_dicts)
        #     ba.is_verb_fromSpacy(token, self.histogram_dicts)

        # for word in nltk.word_tokenize(sentence):
        #     ba.is_verb_fromListNLTK(word, self.verbs, self.histogram_dicts)

        # ba.find_SVOs(doc, self.histogram_dicts)
        pass

    def _current_methods(self, pmid, sentence):

        ba.new_method(
            pmid,
            sentence,
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
        for key in self.HIT_TUPLES.keys():
            try:
                shutil.copyfile(os.path.join(STATS_OUT_DIR, f'hits {key}.txt'), os.path.join(
                    STATS_OUT_DIR, f'hits {key} - BACKUP.txt'))
            except:
                print(
                    f'no backup of {key} made as no previous stats were stored')
            with open(os.path.join(STATS_OUT_DIR, f'hits {key}.txt'), 'w+'):
                pass
                print('backupped and then emptied outputfile:', os.path.join(
                    STATS_OUT_DIR, f'hits {key}.txt'))

    def output_one_pmid_hit(self, hit):
        with open(os.path.join(STATS_OUT_DIR, f'pmid hits.txt'), 'a+') as fh:
            fh.write(f'{hit}\n')

    def _output_hitstrings(self):
        
        all_occurences = Counter()
        for key, tuples in self.HIT_TUPLES.items():
            
            occurrences_small = Counter()
            with open(os.path.join(STATS_OUT_DIR, f'all hits {key}.txt'), 'a+') as fh:
                    for pmid, struct, verb, disfunc_name in tuples:
                        fh.write(
                            f'{pmid};;; {struct};;; {verb};;; {disfunc_name}\n')
                        
                        counter_key = (s.lower() for s in [struct, verb, disfunc_name])
                        occurrences_small[counter_key]+=1
                        all_occurences[counter_key]+=1

            with open(os.path.join(STATS_OUT_DIR, f'hits {key}.txt'), 'a+') as fh:
                for (s,v,f), count in occurrences_small.most_common():
                    fh.write(f'{s}\t{v}\t{f}\t{count}')

        with open(os.path.join(STATS_OUT_DIR, f'hits all-counts.txt'), 'a+') as fh:
            for (s,v,f), count in all_occurences.most_common():
                fh.write(f'{s}\t{v}\t{f}\t{count}')

    def output_stats(self):
        print('outputting stats . . .')

        self._output_histograms()
        self._output_hitstrings()
        self._output_exclude_reasons()

        print('done outputting stats')

    def _output_exclude_reasons(self):
        stats_title_string = f"## STATS EXCLUDE REASONS!"
        headers_string = 'reason, count'
        total_records_string = f"## records analysed for this data: {self.global_stats['records analysed']}"
        total_excluded = sum(val for val in self.excluded_by_counts.values())
        total_excluded_string = f"## records excluded for this data: {total_excluded}"
        print(stats_title_string)
        print(total_records_string)
        with open(os.path.join(STATS_OUT_DIR, f'Exclude reasons counts'), 'w+') as fh:
            fh.write(f"{stats_title_string}\n\n")
            fh.write(f"{total_records_string}\n")
            fh.write(f"{total_excluded_string}\n\n")
            fh.write(f"{headers_string}\n")
            for reason, count in self.excluded_by_counts.items():
                fh.write(f'{reason}, {count}\n')

                with open(os.path.join(STATS_OUT_DIR, f'excluded by {reason} pmids-list'), 'w+') as fh2:
                    for pmid in self.excluded_by_pmids[reason]:
                        fh2.write(f'{pmid}\n')

    def _prep_refactor_list(self, listname, list):
        print(f'writing list about {listname} to file')
        with open(os.path.join(DATA_DIR, f'(refactor)-{listname}'), 'w+') as fh:
            for item in sorted(list):
                fh.write(f'{item}\n')

    def prepare_refactor(self):
        """Since all lists are available here, output them all to files
        to be used after refactoring"""
        self._prep_refactor_list('structures', self.BRAIN_STRUCTURES)
        self._prep_refactor_list('verbs', self.verbs)
        self._prep_refactor_list('functions', self.BRAIN_FUNCTIONS)
        self._prep_refactor_list('disorders', self.DISORDERS)
        self._prep_refactor_list('mesh_include', self.mesh_include_list)
        self._prep_refactor_list('mesh_exclude', self.mesh_exclude_list)
        self._prep_refactor_list('verbs_excluded', self.verbs_excluded)
        self._prep_refactor_list('verbs_included', self.verbs_included)


    def _output_histograms(self):
        for name, counter in self.histogram_dicts.items():
            stats_title_string = f"## STATS FOR {name.upper()} !"
            headers_string = 'item, occurences'
            total_records_string = f"## records analysed for this data: {self.global_stats['records analysed']}"
            print(stats_title_string)
            print(total_records_string)

            with open(os.path.join(STATS_OUT_DIR, f'stats_{name}.txt'), 'w+') as fh:
                fh.write(f"{stats_title_string}\n\n")
                fh.write(f"{total_records_string}\n\n")
                fh.write(f"{headers_string}\n")
                for item, count in counter.most_common():
                    fh.write(f'{item}, {count}\n')


if __name__ == '__main__':
    m = MedlineAnalyser()

    m.prepare_refactor()
    print('reminder: disabled run at bottom of file ')
    # m.run()
