from __future__ import print_function
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import LatentDirichletAllocation
import re, time
from multi_preprocess import *
from collections import defaultdict, Counter

logging.basicConfig(filename='.app.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


# def get_data_and_ner(pmid):
#     biodocs = retrieveBioDocs(str(pmid))
#     data_samples, neslist, total_sentences, sum_tokens = loadBioDoc(biodocs)
#     return data_samples
#
# def buildDict(data_samples):
#     wordsDict = defaultdict(lambda: 0)
#     for doc in data_samples:
#       words = doc.split(" ")
#       for w in words:
#         wordsDict[w] += 1
#     return wordsDict
#
# # filter out dates and websites
# def filter_data(data_samples, wordsDict):
#   total_words = len(wordsDict.keys())
#   top_percent = total_words * .05 #filter top 3& occuring words
#   bottom_percent = total_words * .90 #filter lowest 20% occuring words
#   ordered_dict = Counter(wordsDict)
#   top_words = ordered_dict.most_common(int(top_percent))
#   bottom_words = ordered_dict.most_common()[ int(bottom_percent):  ]
#
#   exclude_top = [ t[0] for t in top_words]
#   exclude_bottom = [b[0] for b in bottom_words]
#   exclude_words = ['http', 'www']
#   all_exclude = exclude_top + exclude_bottom + exclude_words
#
#   filter_data_samples = []
#   for doc in data_samples:
#     temp_doc = []
#     doc_words = doc.split(" ")
#     for word in doc_words:
#       if word not in all_exclude:
#         temp_doc.append(word)
#
#     filter_data_samples.append( " ".join(temp_doc)  )
#   return filter_data_samples




def get_tfidf(data): #data should be a list of strings for the documents
  logging.info("* Preparing to vectorize data ...")
  tfidf_vectorizer = TfidfVectorizer(stop_words='english', ngram_range=(2, 4), norm='l2')
  logging.info("* Fitting data to vector ...")
  tfidf = tfidf_vectorizer.fit_transform(data)
  logging.info("* Successfully fit data to the vector !!! ")
  return tfidf, tfidf_vectorizer


def fit_lda(tfidf, num_topics):
  logging.info("* Initializing Latent Dirichlet Allocation ... ")
  lda = LatentDirichletAllocation(n_topics=num_topics, max_iter=10, learning_method='online', learning_offset=50., random_state=3)
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
  jDict = re.sub('\'', '\"', str(jDict)) #json needs double quotes, not single quotes
  return jDict


def topics_lda(tf_vectorizer, lda, n_top_words):
  #print("\nTopics in LDA model:")
  tf_feature_names = tf_vectorizer.get_feature_names()
  jsonLDA = print_top_words(lda, tf_feature_names, n_top_words)
  return jsonLDA


# data_samples = get_data_and_ner("9108111")
# wordsDict = buildDict(data_samples)
# filter_data_samples = filter_data(data_samples, wordsDict)
# tfidf, tfidf_vectorizer = get_tfidf(filter_data_samples)
# num_topics = 10
# lda = fit_lda(tfidf, num_topics)
# n_top_words = 15
# jsonLDA = topics_lda(tfidf_vectorizer, lda, n_top_words)
