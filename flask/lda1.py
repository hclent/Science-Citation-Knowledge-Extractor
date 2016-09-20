from __future__ import print_function
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import LatentDirichletAllocation
import re, time
from multi_preprocess import *
from nltk import ngrams
from collections import defaultdict

# def get_data_and_ner(pmid):
# 	biodocs = retrieveBioDocs(str(pmid)) #a bunch of strings
# 	data_samples, neslist = loadBioDoc(biodocs)
# 	return data_samples, neslist
#
#
# def makeNgrams(data_samples, n):
#   t1 = time.time()
#   #ngramDict = defaultdict(lambda:0)
#   ngram_list = []
#   for sentence in data_samples:
#     sentence = sentence.split()
#     grams = [sentence[i:i+n] for i in range(int(len(sentence)-n))]
#     for g in grams:
#       space = ' '
#       ngram = str(space.join(g))
#       ngram_list.append(ngram)
#   print("Execute everything: done in %0.3fs." % (time.time() - t1))
#   return ngram_list




def get_tfidf(data): #data should be a list of strings for the documents
  print("* Preparing to vectorize data ...")
  tfidf_vectorizer = TfidfVectorizer(stop_words='english', ngram_range=(1, 2), norm='l2')
  print("* Fitting data to vector ...")
  tfidf = tfidf_vectorizer.fit_transform(data)
  print("* Successfully fit data to the vector !!! ")
  return tfidf, tfidf_vectorizer


def fit_lda(tfidf, num_topics):
  print("* Initializing Latent Dirichlet Allocation ... ")
  lda = LatentDirichletAllocation(n_topics=num_topics, max_iter=25, learning_method='online', learning_offset=50., random_state=1)
  lda.fit(tfidf)
  print("* Successfully fit data to the model!!! ")
  return lda


def print_top_words(model, feature_names, n_top_words):
  jDict = {"name": "flare", "children": []} #initialize dict for json
  for topic_idx, topic in enumerate(model.components_):
    print("Topic #%d:" % topic_idx)
    running_name = 'concept'+str(topic_idx)
    concept_Dict = {"name": running_name, "children": []}
    jDict["children"].append(concept_Dict)
    #print(", ".join([feature_names[i] for i in topic.argsort()[:-n_top_words - 1:-1]]))
    topic_list = ([feature_names[i] for i in topic.argsort()[:-n_top_words - 1:-1]])
    for term in topic_list:
      print(term)
      # if term in ngram_list:
      #   print("TERM IS IN NGRAM_LIST")
      #   for ngram in ngram_list:
      #     if ngram == term:
      #       print(" +1 ")
      # for sentences in data_samples:
      #   sents = sentences.split(" ")
      #   if term in sents:
      #     print("TERM IS IN SENTS (1gram)")
      #     for word in sents:
      #       if word == sents:
      #         print(" + 1 ")
      term_Dict = {"name": term, "size": 700}
      concept_Dict["children"].append(term_Dict)
  jDict = re.sub('\'', '\"', str(jDict)) #json needs double quotes, not single quotes
  return jDict


def topics_lda(tf_vectorizer, lda, n_top_words):
  print("\nTopics in LDA model:")
  tf_feature_names = tf_vectorizer.get_feature_names()
  jsonLDA = print_top_words(lda, tf_feature_names, n_top_words)
  return jsonLDA


# data_samples, nes_list = get_data_and_ner(18269575)
# ngram_list = makeFreqDict(data_samples, 2)
# tfidf, tfidf_vectorizer = get_tfidf(data_samples)
# num_topics = 5
# lda = fit_lda(tfidf, num_topics)
# n_top_words = 5
# jsonLDA = topics_lda(tfidf_vectorizer, lda, n_top_words)
# print(jsonLDA)
