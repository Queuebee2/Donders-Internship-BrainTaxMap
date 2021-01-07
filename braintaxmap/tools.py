import functools
import os
import time
from os import path, stat
from random import random
from time import sleep

FUZZING = True
DATA_DIR = ".." + os.sep + 'data' + os.sep


class DataLoggerHelper():
    threshold = 2

    @staticmethod
    def log_topmesh(fh,  query, stats, top=10):
        fh.write(f'\n-=-=-= TOP 10 MeSH FOR {query} =-=-=-=-\n')
        for term, count in stats['meshterm_occurences'].most_common(top):
            if count > DataLoggerHelper.threshold:
                fh.write(f'{term}, {count}\n')

    @staticmethod
    def log_lastmesh(fh, query, stats, last=10):
        fh.write(f'\n-=-=-= LAST 10 MeSH FOR {query} =-=-=-=-\n')
        for term, count in stats['meshterm_occurences'].most_common()[:-last:-1]:
            if count > DataLoggerHelper.threshold:
                fh.write(f'{term}, {count}\n')

    @staticmethod
    def log_queries_sent(fh, stats):
        queries = ', '.join(stats['queries'].keys())
        fh.write(f"QUERIES SEARCHED:\n{queries}")

    @staticmethod
    def log_allmesh(fh, stats):
        fh.write('\n-=-=-= ALL MeSHterm occurences -=-=-=\n')
        for term, count in stats['allMeshtermCounts'].most_common():
            if count > DataLoggerHelper.threshold:
                fh.write(f'{term}, {count}\n')

    @staticmethod
    def log_pmid_SVOs(fh, query,  stats):
        fh.write(f'\n-=-=-= pmid and SVOs found for {query} =-=-=-=-\n')
        for pmid, svos in stats['svo_pmids'].items():
            fh.write(f"{pmid}, {svos}\n")

    @staticmethod
    def log_all_svos(fh, stats):

        fh.write('\n-=-=-= ALL SVO occurences -=-=-=\n')
        try:
            for svo, count in stats['allSVOs'].most_common():
                if count > DataLoggerHelper.threshold:
                    fh.write(f"{svo[1:-1]}  {count}\n")
        except KeyError:
            fh.write(f'KeyError - No svos?\n')

    @staticmethod
    def log_all_nounphrases(fh, stats):
        fh.write('\n-=-=-= ALL NOUNPHRASE occurences -=-=-=\n')
        for nounphrase, count in stats['nounphrases'].most_common():
            if count > DataLoggerHelper.threshold:
                fh.write(f"{nounphrase}  {count}\n")

    @staticmethod
    def log_all_verbs(fh, stats):
        fh.write('\n-=-=-= ALL VERB occurences -=-=-=\n')
        for verb, count in stats['verbs'].most_common():
            if count > DataLoggerHelper.threshold:
                fh.write(f"{verb}  {count}\n")

    @staticmethod
    def log_all_custom_verbs(fh, stats):
        fh.write('\n-=-=-= ALL (CUSTOM) VERB occurences -=-=-=\n')
        for verb, count in stats['custom_verbs'].most_common():
            if count > DataLoggerHelper.threshold:
                fh.write(f"{verb}  {count}\n")

    @staticmethod
    def log_abstracts(fh, records, meshterms=True,
                      svos=True, meshterm_amnt=0, svo_amt=0):
        """prints a list of selected records (list should already
        have pre-selected records, all will be printed
        """
        for rec in records:
            fh.write(f"\n---------")
            fh.write(f"\nPMID: {rec['PMID']}\n")
            fh.write(f"title: {rec['TI']}\n")
            fh.write(f"ABSTRACT:\n{rec['AB']}\n")

            if meshterms:
                fh.write("\nMeshTerms\n")
                if 'MH' in rec.keys():
                    for i, m in enumerate(rec['MH']):
                        if meshterm_amnt > 0 and i < meshterm_amnt or meshterm_amnt == 0:
                            fh.write(f' - mesh {i}: {m}\n')

            if svos:
                fh.write("\nFound SVOs:\n")
                for i, svo in enumerate(rec['SVOs']):
                    if svo_amt > 0 and i < svo_amt or svo_amt == 0:
                        fh.write(f" - {i+1}. {svo}\n")


def _create_verblist():

    all_words = set()
    with open(DATA_DIR + 'verbs') as fh:
        c = 0
        for l in fh:
            line = l.strip()
            c += 1
            cells = line.split('\t')
            if len(cells) == 6:
                words = cells[1:]
                # print(words)
                for word in words:
                    if '/' in word:
                        a, b = word.replace(" ", '').split('/')
                        all_words.add(a)
                        all_words.add(b)
                        continue
                    elif ('(') in word:
                        word = word.split(' ')[0]
                    all_words.add(word)

    print(len(all_words))
    print(c)
    with open(DATA_DIR + '1000-verbs-set.txt', 'w+') as fh:
        fh.writelines([f'{verb}\n' for verb in all_words])


def load_verbs(filepath=DATA_DIR + '1000-verbs-set.txt'):
    with open(filepath, 'r') as fh:
        verbs = set()
        for line in fh:
            verbs.add(line.strip())
    print(f'successfuly loaded {len(verbs)} custom verbs')
    return verbs


def fuzz(amt=False, fuzzing=FUZZING):
    """Inspired by Raymond Hettingers python talks
    https://github.com/rhettinger
    """
    if fuzzing:
        if not amt:
            sleep(random())
        else:
            sleep(amt)


def timethisfunc(func):
    """wrapper adapter from 
    https://realpython.com/lessons/timing-functions-decorators/
    and 
    https://medium.com/pythonhive/python-decorator-to-measure-the-execution-time-of-methods-fa04cb6bb36d
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ts = time.perf_counter()
        result = func(*args, **kwargs)
        te = time.perf_counter()
        duration = time.strftime('%Mm%Ss', time.gmtime(te-ts))
        print(f"{func.__name__!r} ran for {duration}")
        return result
    return wrapper


def centr_print(title='', motif='- ', amt=60):
    """center some text with custom side motif"""
    title = title + ' ' if len(title) % 2 == 0 else title
    charsleft = abs(amt - len(title))
    side = (charsleft // 2 // len(motif))
    print(f"{(side * motif)}{motif[0]} {title} {side * motif}")


def tester(TESTS):
    """

    @param TESTS:  dict of {'program':(mode, funcname)}
    """
    if len(TESTS) == 0:
        print('No programs to run!')
        return
    centr_print(f'Testing stuff!')

    for testkey, (testmode, testfunc) in TESTS.items():
        if testmode:
            centr_print(f'testing {testkey}')
            testfunc()
            centr_print(f'Done testing {testkey}')
        else:
            print(f'testing is DISABLED for {testkey}')

    centr_print(f'All testing done')


def programrunner(PROGRAMS):
    """

    @param PROGRAMS:  dict of {'program':(mode, funcname)}
    """
    if len(PROGRAMS) == 0:
        print('No programs to run!')
        return
    centr_print(f'Running programs')

    for progkey, (mode, func) in PROGRAMS.items():
        if mode:
            centr_print(f'running {progkey}')
            func()
            centr_print(f'Done testing {progkey}')
        else:
            print(f'run is DISABLED for {progkey}')

    centr_print(f'Done running programs')


if __name__ == '__main__':
    print(f'running from braintaxmap.tools.py as __main__')
