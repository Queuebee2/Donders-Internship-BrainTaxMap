from neo4j_article_inserts_automated import harvest_articles, grab_keywords
from nltk.corpus import wordnet, stopwords
from nltk.tokenize import word_tokenize
from querymachine import QueryMachine

# working with this as reference
# https://colab.research.google.com/drive/1rnTjTjZKSP-pcRPzCu5mBBZAHZ4DG8yJ#scrollTo=rO98sY97ZzuR


import re
punctionation = re.compile(r'[-.?!,:;()|0-9]')

def fix(rec):
    if 'MH' not in rec.keys():
        rec['MH'] = []
    return rec

def bigprint(rows=1, *objs):
    for r in range(rows):
        print(45*'=')
    for t in objs:
        print(t)
    for r in range(rows):
        print(45*'=')

def test_similarities():
    behaviors, structures = grab_keywords()
    articles = harvest_articles(100, True, 100)
    print([x for x in behaviors][:5])
    print([x for x in structures][5:])
    filter_words = stopwords.words('english')
    for (word, labels), rec in articles:
        print(word, labels, rec)

        # grab abstract and mesh terms
        abstract = rec['AB']
        mesh_terms = fix(rec)['MH']

        # tokenize the abstract
        ab_tokens = word_tokenize(abstract)

        # filter tokens -> remove stopwords, remove punctuations
        post_filter = []
        for words in ab_tokens:
          word = punctionation.sub("", words)
          if len(word)>0 and word not in filter_words:
            post_filter.append(word)


        print(post_filter[:5])


        # determine similarity... HOW
        # https://www.nltk.org/howto/wordnet.html ?


def is_relevant_article(tokens, lists):
    """ based on a list of tokens and a list of lists, if any token
    appears in any of the lists, return True, else, return False
    """
    return any([token in l for l in lists for token in tokens])


veto_words = ['cortex', 'mice', 'mouse']

def is_relevant_2(abstract):
    for word in veto_words:
        if word in abstract:
            return True
    return False

if __name__ == '__main__':
    behaviors, structures = grab_keywords()
    keywords = ['barrel cortex', 'barrel field']
    filter_words = stopwords.words('english')
    q = QueryMachine
    no_synsets = []

    for search_word in keywords:
        records = q.queryPubMed(search_word, 10000)

        for rec in records:
            print(rec['TI'])
            if 'AB' not in rec.keys() or not rec['AB']:
                print('no abstract')
                continue
            abstract = rec['AB']
            mesh_terms = fix(rec)['MH']

            ab_tokens = word_tokenize(abstract)

            # filter tokens -> remove stopwords, remove punctuations
            post_filter = []
            for words in ab_tokens:
                word = punctionation.sub("", words)
                if len(word) > 0 and word not in filter_words:
                    post_filter.append(word)

            # skip if no terms in the abstract occur in our lists
            if not is_relevant_article(post_filter, [behaviors, structures]):
                print('not relevant:',end= ' ')
                print(post_filter)
                # second check...
                if not is_relevant_2(abstract):
                    print('really not relevant')
                    continue

            # check similarities between behaviours, structures and found new words
            for word in post_filter:


                # take only new words
                if word not in behaviors and word not in structures:
                    # print('scores:')
                    word_synsets = wordnet.synsets(word)
                    if not word_synsets:
                        print(word,'has no synsets')
                        continue
                    word_syn = word_synsets[0]

                    scores = {'behave_path':0,
                              'struct_path':0,
                              'behave_lch':0,
                              'struct_lch':0,
                              }

                    for x in behaviors:
                        syn2 = wordnet.synsets(x)[0]
                        path = word_syn.path_similarity(syn2)
                        lch = word_syn.lch_similarity(syn2)
                        scores['behave_path'] += path
                        scores['behave_lch'] += lch

                        if path > 0.8 or lch > 0.8:
                            print('high sim:', path, lch, syn2, word_syn)

                    for x in structures:
                        syn2 = wordnet.synsets(x)[0]
                        path = word_syn.path_similarity(syn2)
                        lch = word_syn.lch_similarity(syn2)
                        scores['struct_path'] += path
                        scores['struct_lch']  += lch
                        if path > 0.8 or lch > 0.8:
                            print('high sim:', path, lch, syn2, word_syn)

                    # print(scores)
                    # We need to determine a clear set of rules.
                    # how do we score a word, what does the score mean
                    # Do we look for SEMANTIC similarity?

                    # will we make/train a classifier? that picks the right words ?
                    # https: // www.tensorflow.org / tutorials / text / word_embeddings

                    # https://spacy.io/usage/vectors-similarity




