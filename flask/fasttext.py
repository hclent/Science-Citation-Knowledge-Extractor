from multi_preprocess import retrieveBioDocs
from processors import *
import sys, os.path, time
import numpy as np
from collections import defaultdict, Counter
import gensim, os, codecs
from gensim.models import Word2Vec
from sklearn.cluster import KMeans
import math, random
from kmeans1 import do_kemeans, get_hashing


'''
Step 1:
    Get lemmas and tags from biodocs
    [[doc1 words], [doc2 words]] & [[doc1 tags],[doc2 tags]]
Step 2:
    NP tokenization
    [[doc1 tokens], [doc2 tokens]]
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


def flatten(listOfLists):
    return list(chain.from_iterable(listOfLists))

#Input: biodoc
#Output, list of words, and corresponding list of tags
def get_words_tags(pmid_list):
    words = []
    tags = []
    for pmid in pmid_list:
        biodocs = retrieveBioDocs(pmid)
        #print(len(biodocs))
        for doc in biodocs:
            doc_words = []
            doc_tags = []
            json_file = doc["jsonpath"]
            with open(json_file) as jf:
                data = Document.load_from_JSON(json.load(jf))
                num_sentences = data.size
                for i in range(0, num_sentences):
                    s = data.sentences[i]
                    s_words = s.lemmas
                    s_tags = s.tags
                    doc_words.append(s_words)
                    doc_tags.append(s_tags)
            doc_words = flatten(doc_words)
            doc_tags = flatten(doc_tags)
            words.append(doc_words)
            tags.append(doc_tags)
    flat_words = flatten(words)
    flat_tags = flatten(tags)
    return flat_words, flat_tags

#Input: list of words and corresponding list of tags
#Output: list of noun phrases joined as a string
def transform_text(words, tags):
    transformed_tokens = []
    np_stack = []
    join_char = " " #instead of ("_")
    for i, w in enumerate(words):
        if tags[i].startswith("NN"):
            np_stack.append(w)
        else:
            if len(np_stack) > 0:
                np = join_char.join(w.lower() for w in np_stack)
                transformed_tokens.append(np)
                # reset stack
                np_stack = []
            # append current token
            transformed_tokens.append(w)
    return transformed_tokens

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



