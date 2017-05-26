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
            with open(doc) as jf:
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






#pmid_list = retrieveAllPmids()
# pmid_list = ['18952863','18269575']
# words, tags = get_words_tags(pmid_list) #list of words/tags per doc
# transformed_sentence = transform_text(words, tags)
# npDict = chooseTopNPs(transformed_sentence)
# top_nps= list(npDict.most_common(100))
# # ### top_nps200 = list(npDict.most_common(150))
# # ### top_nps = [item for item in top_nps200  if item not in top_nps100]
# # print(top_nps)
# #
# model = load_model('/home/hclent/tmp/fastText/16kmodel.vec')
# print(model.syn0.shape)
# matrix = getNPvecs(top_nps, model)
# kmeans = KMeans(n_clusters=10, random_state=2).fit(matrix)
# #
# #
# res = list(zip(kmeans.labels_, top_nps))
# print(res)
# #
# t1 = []
# t2 = []
# t3 = []
# t4 = []
# t5 = []
# t6 = []
# t7 = []
# t8 = []
# t9 = []
# t10 = []
# for t in res:
#     if t[0] == 0:
#         t1.append(t[1])
#     if t[0] == 1:
#         t2.append(t[1])
#     if t[0] == 2:
#         t3.append(t[1])
#     if t[0] == 3:
#         t4.append(t[1])
#     if t[0] == 4:
#         t5.append(t[1])
#     if t[0] == 5:
#         t6.append(t[1])
#     if t[0] == 6:
#         t7.append(t[1])
#     if t[0] == 7:
#         t8.append(t[1])
#     if t[0] == 8:
#         t9.append(t[1])
#     if t[0] == 9:
#         t10.append(t[1])
# print("t1: ")
# print(t1)
# print("t2: ")
# print(t2)
# print("t3: ")
# print(t3)
# print("t4: ")
# print(t4)
# print("t5: ")
# print(t5)
# print("t6: ")
# print(t6)
# print("t7: ")
# print(t7)
# print("t8: ")
# print(t8)
# print("t9: ")
# print(t9)
# print("t10: ")
# print(t10)
# print("results: ")
# results = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10]
# print(results)





#################  graveyard ###########################
# def makeFeatures(words, tags, top_nps, model):
#     feature_vecs = []
#     size = model.syn0.shape
#     num_features = size[1]
#     print(top_nps)
#     i = 0
#     for doc in words:
#         doc_vec = []
#         transformed = transform_text(doc, tags[i])
#         #print(transformed)
#         for noun in top_nps: #top_ns is a counter
#             #print(noun)
#
#             if noun[0] in transformed:
#                 np_vec = model[noun[0]]
#                 #print(np_vec)
#                 doc_vec.append(np_vec)
#             else:
#                 v_list = [0.0] * int(num_features)
#                 empty_vec = np.array(v_list)
#                 #print(empty_vec)
#                 doc_vec.append(empty_vec)
#                 #append a vec of length num_features
#         doc_vec = np.array(flatten(doc_vec))
#         feature_vecs.append(doc_vec)
#
#         i += 1
#
#     print(len(feature_vecs))
#     print(type(feature_vecs[0]))
#     print(feature_vecs[0])
#     print(len(feature_vecs[0]))
#     return feature_vecs


# fv = makeFeatures(words, tags, top_nps, model)
# print(fv)
#
#
# hX, hasher = get_hashing(fv) #returns hX and hasher
# clusters = do_kemeans(hX, 10) #returns clusters
#
# rank = zip(list(top_nps, clusters))
# println(rank)

#### For alternative tokenization  #### noun_phrase_words
#
# print("making transformed s")
# transformed_s = (' '.join(transformed_sentence))
# print(transformed_s[:1000])
# suffix = "/home/hclent/data/nns/5k"
# completeName = os.path.join(suffix, ('5210_nns.txt'))  #pmcid.txt #save to suffix path
# sys.stdout = open(completeName, "w")
# print(transformed_s)

##model = load_model('/home/hclent/tmp/fastText/ftmodel.vec') #fasttext trained on 5k
##model = load_model("/home/hclent/data/nns/5k/5210_nns") #wtv trained on 5k