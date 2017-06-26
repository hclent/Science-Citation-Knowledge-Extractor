from __future__ import print_function
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import LatentDirichletAllocation
import json, logging


logging.basicConfig(filename='.app.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


def get_tfidf(data): #data should be a list of strings for the documents (lemmas_for_lda)
  list_data_strings = [' '.join(map(str, d)) for d in data] #map to string. strings are necessary for the TFIDF
  logging.info("* Preparing to vectorize data ...")
  tfidf_vectorizer = TfidfVectorizer(stop_words='english', ngram_range=(2, 4), norm='l2')
  logging.info("* Fitting data to vector ...")
  tfidf = tfidf_vectorizer.fit_transform(list_data_strings)
  logging.info("* Successfully fit data to the vector !!! ")
  return tfidf, tfidf_vectorizer


#n_jobs for speed :)
#TODO: cannot set n_jobs above 1. Warning: Multiprocessing backend parallel loops cannot be nested below threads
def fit_lda(tfidf, num_topics):
  logging.info("* Initializing Latent Dirichlet Allocation ... ")
  lda = LatentDirichletAllocation(n_topics=num_topics, max_iter=10, learning_method='online', learning_offset=50., random_state=3, n_jobs=5)
  lda.fit(tfidf)
  logging.info("* Successfully fit data to the model!!! ")
  return lda


def print_top_words(model, feature_names, n_top_words):
  jDict = {"name": "flare", "children": []} #initialize dict for json
  for topic_idx, topic in enumerate(model.components_):
    #print("Topic #%d:" % topic_idx)
    running_name = 'concept'+str(topic_idx)
    concept_Dict = {"name": running_name, "children": []}
    jDict["children"].append(concept_Dict)
    #print(", ".join([feature_names[i] for i in topic.argsort()[:-n_top_words - 1:-1]]))
    topic_list = ([feature_names[i] for i in topic.argsort()[:-n_top_words - 1:-1]])
    for term in topic_list:
      # print(term)
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

  jsonDict = json.dumps(jDict)

  return jsonDict


def topics_lda(tf_vectorizer, lda, n_top_words):
  #print("\nTopics in LDA model:")
  tf_feature_names = tf_vectorizer.get_feature_names()
  jsonLDA = print_top_words(lda, tf_feature_names, n_top_words)
  return jsonLDA
