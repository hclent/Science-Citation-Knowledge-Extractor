from __future__ import print_function
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.decomposition import LatentDirichletAllocation
import re


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
      term_Dict = {"name": term, "size": 700}
      concept_Dict["children"].append(term_Dict)
  jDict = re.sub('\'', '\"', str(jDict)) #json needs double quotes, not single quotes
  return jDict


def topics_lda(tf_vectorizer, lda, n_top_words):
  print("\nTopics in LDA model:")
  tf_feature_names = tf_vectorizer.get_feature_names()
  jsonLDA = print_top_words(lda, tf_feature_names, n_top_words)
  return jsonLDA


# tfidf, tfidf_vectorizer = get_tfidf(data_samples)
# lda = fit_lda(tfidf, num_topics)
# jsonLDA = topics_lda(tfidf_vectorizer, lda, n_top_words)
# print(jsonDict)
