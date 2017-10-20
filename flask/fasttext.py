from processors import *
import numpy as np
import pickle, os, re
from collections import defaultdict, Counter, Iterable
from gensim.models import Word2Vec
from configapp import app


'''
Step 1:
    Get lemmas and tags from cache
Step 2:
    NP tokenization
    [[doc1 noun tokens], [doc2 noun tokens]]
Step 3:
    flatmap tokens for training w2v
    [all tokens ever]
    train w2v on all tokens
Step 4:
    Try these:
    Get the top X% NPs,
    Rank NPs with TF-IDF, Or just get the top x NPs
Step 5:
    cluster
'''


vis_stop_words = ['university', '%', 'table', 'figure', '\\u', '\\\\', '\\', 'author', 'publication', 'appendix',
                      'table', 'author', 'skip', 'main', '.', 'title', 'u2009', 'publisher', 'article',
                      'www.plantphysiol.org', 'copyright', 'san diego', 'california',  '. . . .', '. .', ',', '.....',
                      "\"", "1", ";", "3", '.', ' . ',
                      "one", "also", "=", "2", "4" "number", 'j.', 'm.', 's.', 'many', 'b', '6', '10', 'however',
                      'well', 'c', 'p.', '*', "'s", ':', "'", '0', '4', '-', 'three', 'may', 'non', 'could',
                      'would', 'two', 'one', '.',
                      'e.g.', 'doi', 'case', 'follow', 'describe', 'name', 'see', 'among', 'single', 'several',
                      'run', 'additional','number', 'show', 'include', 'use', 'multiple', 'important', 'individual', 'like',
                      'exist', 'related', 'gateway', 'control', 'suggest', 'high', 'allow', 'first', 'list', 'define', 'set',
                      'new', 'different', 'thus', 'small', 'year', 'due', 'i.e.', '...', 'low', 'per', 'big', 'via',
                      '. . . . .', 'full', 'another', 'second', '100', 'google', 'u2003', 'character', 'state', 'ieee',
                      'july', 'vol', 'username', 'password', 'email', 'address', 'e01797', 'swetnam', '1 0/27', '1 0/21', 'sen.',
                      'license', 'biorxiv', 'preprint', 'textgoogle', 'crossrefpubmedgoogle', 'john', 'wiley', 'peerj', 'apr.', 'view', 'scopus',
                      'crossrefpubmedweb', 'sciencegoogle', 'scholar', '|', 'journal', 'apus', 'aug', 'publ', 'jan.', 'jan', 'sep',
                      'springer', 'rights', 'elsevier', 'b.v.', 'appl19', 'biol', 'oxford', 'press', 'viewerdownload', 'powerpoint',
                      'librarycopyright', 'letter', 'interest', 'tip', 'note', 'april', 'volume', 'isbn', "{", "}" "[", "]", 'n', 'user' ]



def flatten(listOfLists):
    return list(chain.from_iterable(listOfLists))


#Input: query
#Output, list of words, and corresponding list of tags
def get_words_tags(query):
    pmid_list = query.split('+')
    first_pmid = pmid_list[0]


    #path_to_lemma_samples = '/home/hclent/data/pmcids/' + str(first_pmid[0:3]) + '/' + str(
    #    first_pmid[3:6]) + '/lemma_samples_' + str(query) + ".pickle"
    name = str(first_pmid[0:3]) + '/' + str(first_pmid[3:6]) + '/lemma_samples_' + str(query) + ".pickle"

    path_to_lemma_samples = os.path.join((app.config['PATH_TO_CACHE']), name)
    with open(path_to_lemma_samples, "rb") as file:
        lemma_samples = pickle.load(file)

    words = [l[1] for l in lemma_samples]

    flat_words = flatten(words)

    tags = [l[2] for l in lemma_samples]
    flat_tags = flatten(tags)

    keep_words = []
    keep_tags = []

    for c in list(zip(flat_words, flat_tags)):
        w = c[0]
        tag = c[1]
        if w not in vis_stop_words and not re.match('^(http|u\d{4})', w):
            keep_words.append(w)
            keep_tags.append(tag)

    word_length = len(keep_words)
    tag_length = len(keep_tags)

    # It SHOULD NOT happen that the tags and words are not the same length...
    # But JUST IN CASE.... force them to have the same length
    if word_length != tag_length:
        #print("they aren't the same :(( ")
        if word_length > tag_length:
            #shorten word_length
            keep_words = keep_words [:tag_length]
        elif tag_length > word_length:
            #shorten tag_length
            keep_tags = keep_tags[:word_length]

    return keep_words, keep_tags




#Input: list of words and corresponding list of tags
#Output: list of noun phrases joined as a string
def transform_text(words, tags):
    transformed_tokens = []
    np_stack = []
    join_char = " " #instead of ("_")
    for i, w in enumerate(words):
        if tags[i].startswith("NN"):
            #print(w)
            #print(i)
            np_stack.append(w)
        else:
            if len(np_stack) > 0:
                np = join_char.join(w.lower() for w in np_stack)
                transformed_tokens.append(np)
                # reset stack
                np_stack = []
            # append current token
            transformed_tokens.append(w)
    keep_transformed_tokens = [w for w in transformed_tokens if len(w.split()) < 4]
    return keep_transformed_tokens



#make dict of NPs with their frequency counts
def chooseTopNPs(transformed_tokens):
    npDict = Counter(defaultdict(lambda: 0))
    for t in transformed_tokens:
        if " " in t: #for previous tokenization strategy, if "_"
            npDict[t] += 1
    return npDict


#Load fasttext model
def load_model(file):
    if ".vec" not in file:
        model = Word2Vec.load(str(file))
    if ".vec" in file:
        model = Word2Vec.load_word2vec_format(file, binary=False)  # C binary format
    model.init_sims(replace=True)
    return model


# make list of embeddings to cluster [[embedding 1], [embedding 2], ... ]
def getNPvecs(top_nps, model):
    all_vecs = []
    for nps in top_nps:
        try:
            nouns = nps[0].split(' ')
            #print(nouns)
            denomintor = len(nouns)
            zeroes = [0.0] * 100
            sum = np.asarray(zeroes).T #init empty column vector (100,)
            for n in nouns:
                try:
                    np_vec = model[n] #(100,)
                    sum += np.add(sum, np_vec)

                except Exception as oov:
                    pass
            # Average over the embeddings (e.g.  sum 2 NP embeddings and divide by 2)
            ### Problem: if one word is oov in the noun phrase, it will still average the embeddings :/
            avg_vec = np.divide(sum, denomintor)
            vec = avg_vec.tolist()
            all_vecs.append(vec)
            #print('-----------------------')
        except Exception as e:
             pass
    matrix = np.array(all_vecs)
    return matrix



