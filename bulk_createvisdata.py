import os
from braintaxmap.tools import (read_icd11_flat, dsm5parse, read_icd11_lvls)
import logging
from bulkconstants import DATA_DIR, DISORDER_RESULT_FILENAMES
from collections import defaultdict, Counter
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s: %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)
import tqdm


def compare_disorders():
    """check whether dsm5 
    are subset of icd11-, otherwise, we have to remove all
    hits from dsm5 that dont occur in icd11 . 
    if we just return the pruned icd11 here thats got to be fine.
    """
    logger.info('comparing disorders to dsm5')
    icd11 = read_icd11_flat()
    dsm5 = dsm5parse()

    different = set()
    for term in dsm5:
        if term.lower() not in icd11:
            logger.info(f'different: {term}')
            different.add(term)

    logger.info(f'different: {len(different)}')
    logger.info(f'DSM5: {len(dsm5)} ICD: {len(icd11)}')
    return different


def get_dsm5_terms_to_filter_out():
    return compare_disorders()


def find_unique_dsm_hits():
    RESULT = "\n"+40*"-"+\
    """
    (already ran)
    PRECALCULATED INFORMATION (find_unique_dsm_hits()):
    there are 2262 results where a dsm5-only term is used that
    is not included in the icd11. The other 188786 results are fine
    
    UPDATE for 'all hits disorders.txt specifically: 793 problems and 5340 good \n""" + 40*"-"

    logger.info(RESULT)
    return

    not_in_icd11 = compare_disorders()
    logger.info('Got disorders that are unique to dsm5')

    problematic = 0
    noproblem = 0
    for fn in DISORDER_RESULT_FILENAMES:

        path = os.path.join(DATA_DIR, 'stats', fn)

        with open(path) as fh:
            logger.info(f'checking output files for {path} for dsm5 terms')
            for line in fh:
                for word in not_in_icd11:
                    if word.lower() in line.lower():
                        logger.info(line.strip())
                        problematic += 1
                        break
                else:
                    noproblem += 1

    logger.info(f'Found {problematic} problems and {noproblem} good')


def _icd_flat_reader(path=os.path.join(DATA_DIR, 'visualisation_help',
                                       'icd11-used-terms.txt')):
    with open(path) as fh:
        for line in fh:
            yield line.strip().lower()


def _line_level(line):
    l = 0
    while len(line) > 0 and not line[0].isalnum():
        line = line[2:]
        l += 1
    return l


def _strip_line(line):

    while len(line) > 0 and not line[0].isalnum():
        line = line[2:]
    return line


def read_icd_flat_to_filtered_hierarchy():
    """Read the icd-11-used-terms.txt file and parse it,creating a dictionary where
    all sub-terms for each level in the hierarchy are part of a set of terms
    for example:
    
    Any terms that only occur in dsm-5 are stripped out

    LEVEL 0: Mental, behavioural or neurodevelopmental disorders 
    LEVEL 1: - Neurodevelopmental disorders
    LEVEL 2: - - Disorders of intellectual development
    LEVEL 3: - - - Disorder of intellectual development, mild
    LEVEL 2: - - Developmental speech or language disorders
    LEVEL 3: - - - Developmental speech sound disorder

    will yield
    result = {
        0:{'Mental, behavioural or neurodevelopmental disorders':{
                'Neurodevelopmental disorders',
                'Disorders of intellectual development',....
            }
        },
        1:{'Neurodevelopmental disorders':{
                'Disorders of intellectual development',
                'Disorder of intellectual development, mild',....
            },
        2:{ ..}
    }
    """
    dicionaries = defaultdict(dict)
    filter_out = get_dsm5_terms_to_filter_out()

    skipped_by_filter = 0
    for level in range(5):
        current = None
        leveldict = defaultdict(set)

        for line in _icd_flat_reader():
            if _strip_line(line) in filter_out:
                skipped_by_filter += 1
                continue
            sublvl = _line_level(line)
            if sublvl == level:
                major = _strip_line(line)
            elif sublvl > level:
                leveldict[major].add(_strip_line(line))

        dicionaries[level] = leveldict

    for level, maindict in dicionaries.items():
        for subdict_name, subdict in maindict.items():
            logger.info(
                F"{len(subdict): 4d} ITEMS IN LEVEL {level}: {subdict_name} ")
    logger.info(
        f"Skipped {skipped_by_filter} terms because of no-dsm5-filter (0=good)"
    )
    logger.info('Returning hierarchy levels for icd-11 file as dict')
    return dicionaries


def _write_treemap_data(data):
    logger.info("Writing treemap data to files")
    for level, main_counter in data.items():
        logger.info(f'Creating data for treemap for level {level}')
        with open(
                os.path.join(DATA_DIR, 'used-for-internship',
                             f'disorder-counts-lvl{level}-treemap.csv'),
                'w+') as fh:
            fh.write(f'name\tcount\n')
            for item, counts in main_counter.most_common():
                fh.write(f'{item}\t{counts}\n')

    logger.info('Treemaps done. I think.')


def create_treemap_data(disorder_level_hierarchy):
    hierarchy_counts = defaultdict(Counter)
    filter_out = get_dsm5_terms_to_filter_out()

    path = os.path.join(DATA_DIR, 'stats', 'all hits disorders.txt')
    with open(path, 'r') as fh:
        c = 0
        for line in fh:
            pmid, struct, verb, disorder = line.strip().split(';;; ')
            for level, sublists in disorder_level_hierarchy.items():
                for sublistname, lst in sublists.items():
                    c += 1
                    if disorder in lst:
                        hierarchy_counts[level][sublistname] += 1

    logger.info(f'Created treemap data with {c} iterations, I think.')
    log_doubledict(hierarchy_counts)
    _write_treemap_data(hierarchy_counts)
    _write_treemap_hierarchy_again(disorder_level_hierarchy)
    return hierarchy_counts


def log_doubledict(dd):
    logger.info("Logging double dict stuff")

    totals = []
    for km, d in dd.items():
        total = 0
        for ks, x in d.items():
            if isinstance(x, int):
                logger.info(f'{km} {ks} {x}')
                total += x
        totals += [(km, total)]

    logger.info("Also, if these were counts:")
    for km, total in totals:
        logger.info(f" the total for {km} is {total}")
    logger.info(f"Total all: {sum([tot for k,tot in totals])}")


def write_hits_to_tsv():
    logger.info("Writing hits to separate tsv files per combination (struct)x(verb), (disorder)x(verb), (disorder)x(struct)")

    # these are dsm-5 terms that have to be filtered out again
    # since dsm-5 and icd-11 were combined
    filter_out = get_dsm5_terms_to_filter_out()
    pathin = os.path.join(DATA_DIR, 'stats', 'all hits disorders.txt')

    disorder_struct_path = os.path.join(DATA_DIR, 'used-for-internship',
                                        'disorders-structure-dsmfiltered.tsv')
    disorder_verb_path = os.path.join(DATA_DIR, 'used-for-internship',
                                      'disorders-verb-dsmfiltered.tsv')
    struct_verb_path = os.path.join(DATA_DIR, 'used-for-internship',
                                    'struct-verb-dsmfiltered.tsv')

    disorder_struct_counts = Counter()
    struct_distorder_counts = Counter()
    disorder_verb_counts = Counter()
    struct_verb_counts = Counter()
    with open(pathin) as fh:
        with open(disorder_struct_path, 'w+') as fh_out_dis_struct, \
             open(disorder_verb_path,'w+') as fh_out_dis_verb,      \
             open(struct_verb_path,'w+') as fh_out_struct_dis:

            fh_out_dis_struct.write("pmid\tdisorder\tstructure\n")
            fh_out_dis_verb.write("pmid\tdisorder\tverb\n")
            fh_out_struct_dis.write("pmid\tstructure\tdisorder\n")
            for line in fh:
                pmid, struct, verb, disorder = line.lower().strip().split(
                    ';;; ')

                if not disorder in filter_out:
                    fh_out_dis_struct.write(f'{pmid}\t{disorder}\t{struct}\n')
                    fh_out_dis_verb.write(f'{pmid}\t{disorder}\t{verb}\n')
                    fh_out_struct_dis.write(f'{pmid}\t{struct}\t{verb}\n')

                    disorder_struct_counts[f'{disorder};{struct}']+=1
                    struct_distorder_counts[f'{struct};{disorder}']+=1
                    disorder_verb_counts[f'{disorder};{verb}']+=1
                    struct_verb_counts[f'{struct};{verb}']+=1
    
    logger.info('writings counter dictionaries to files (for the so maniestthth time)')
    for counterDictName, theCounter in [
        ('dis-struct', disorder_struct_counts),
        ('struct-dis',struct_distorder_counts),
        ('dis-verb',disorder_verb_counts),
        ('struct-verb',struct_verb_counts)
    ]:
        path=os.path.join(DATA_DIR, 'used-for-internship', f'{counterDictName}.txt')
        logger.info(f'writing {counterDictName} top 10 counts (if counts > 1)')
        with open(path, 'w+') as fout:
            for key, counts in theCounter.most_common():
                main, sub = key.split(';')
                line = f'{main}\t{sub}\t{counts}'
                logger.info(line)
                fout.write(f'{line}\n')
            

def create_disorder_struct_count_map_files():
    filter_out = get_dsm5_terms_to_filter_out()
    pathin = os.path.join(DATA_DIR, 'used-for-internship',
                          'disorders-dsmfiltered.tsv')

    mapping = defaultdict(Counter)
    with open(pathin) as fh:
        _headers = next(fh)
        for line in fh:
            pmid, struct, verb, disorder = line.strip().split('\t')
            mapping[struct][disorder] += 1

    for struct, disordercounter in mapping.items():
        with open('struct-disorder-count {struct}.tsv'):
            pass


def store_hits_in_neo4j():
    print("ALREADY RAN THIS FUNCTION. (struct_node, pmid_node)()")
    logger.info("ARLEADY RAN THIS FUNCTION (struct_node, pmid_node()")
    return
    from py2neo import Graph, Node, Relationship
    from braintaxmap.config import neo4j_db_creds, neo4j_URL

    filter_out = get_dsm5_terms_to_filter_out()
    pathin = os.path.join(DATA_DIR, 'stats', 'all hits disorders.txt')

    graph = Graph(neo4j_URL, auth=neo4j_db_creds)
    logger.info("connected to neo4j")
    with open(pathin) as fh:
        for i, line in enumerate(fh):

            if i % 100 == 0:
                logger.info(f'inserting node {i}')

            pmid, struct, verb, disorder = line.lower().strip().split(';;; ')
            


            if not disorder in filter_out:
                VERB = Relationship.type(verb.upper())
                APPEARED_IN = Relationship.type("APPEARED_IN")

                struct_node = Node('structure', name=struct)
                struct_node.__primarylabel__= 'structure'
                struct_node.__primarykey__ = 'name'

                disorder_node = Node('disorder', name=disorder)
                disorder_node.__primarylabel__= 'disorder'
                disorder_node.__primarykey__ = 'name'

                pmid_node = Node('PMID', name=pmid)
                pmid_node.__primarylabel__= 'PMID'
                pmid_node.__primarykey__ = 'name'
                

                graph.merge(VERB(struct_node, disorder_node))
                graph.merge(APPEARED_IN(struct_node, pmid_node)|APPEARED_IN(disorder_node, pmid_node))

    logger.info("Done inserting nodes into neo4j")
def main():
    # done
    # find_unique_dsm_hits()

    # done
    # data = read_icd_flat_to_filtered_hierarchy()
    # create_treemap_data(data)

    # in progress
    # write_hits_to_tsv()
    # in progress
    # create data: struct:counter[disorder]:count


    # store_hits_in_neo4j()

if __name__ == "__main__":
    main()