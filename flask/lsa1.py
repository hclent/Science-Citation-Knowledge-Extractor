from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn import metrics
from sklearn.decomposition import TruncatedSVD
import pickle, re
from multi_preprocess import *
from collections import defaultdict, Counter

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
#   # print(len(filter_data_samples))
#   # print(data_samples[0])
#   # print("###############################################")
#   # print(filter_data_samples[0])
#   return filter_data_samples


#Input: a list of strings for the document
#Output: TFIDF matrix X, and TfidfVectorizer function
#Output: TFIDF matrix X is a sparse matrix
def get_tfidf(data): 
  tfidf_vectorizer = TfidfVectorizer(stop_words='english', ngram_range=(2, 3), norm='l2')
  tfidf = tfidf_vectorizer.fit_transform(data) 
  return tfidf, tfidf_vectorizer


#Input: TFIDF matrix (X) and TfidfVectorizer
#Output: json dict ready for D3 visualization
#Output: this dict has Latent Semantic Analysis topics (concepts) and topic words
#To Do: Normalization with make_pipeline??? But make_pipeline has no components_
def do_LSA(X, vectorizer, k): 
  lsa = TruncatedSVD(n_components=k, n_iter=100)
  lsa_results = lsa.fit(X)
  terms = vectorizer.get_feature_names()
  jDict = {"name": "flare", "children": []} #initialize dict for json
  for i, comp in enumerate(lsa.components_):
    termsInComp = zip(terms, comp)
    sortedTerms = sorted(termsInComp, key=lambda x: x[1], reverse=True)[:7] #[: number of words desired]
    running_name = 'concept'+str(i)
    concept_Dict = {"name": running_name, "children": []}
    jDict["children"].append(concept_Dict)
    for term in sortedTerms:
      #print(term[0])
      term_Dict = {"name": term[0], "size": 700}
      concept_Dict["children"].append(term_Dict)
    #print("#########################")
  jsonDict = re.sub('\'', '\"', str(jDict)) #json needs double quotes, not single quotes
  #print(jsonDict)
  return jsonDict



# data_samples = get_data_and_ner("18952863")
# wordsDict = buildDict(data_samples)
# filter_data_samples = filter_data(data_samples, wordsDict)
# matrix, vectorizer = get_tfidf(filter_data_samples)
# dictionary = do_LSA(matrix, vectorizer, 10)
