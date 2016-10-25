import naive_makeVecs as makeVecs
import pickle
from database_management import db_citation_urls


#Load from pickled data_samples instead of filename
def loadFromDataSamples(data_samples):
    vecs_list = []

    for document in data_samples:
        vectorCounter = makeVecs.text2vec(document)
        vecs_list.append(vectorCounter)
    return vecs_list


#Load txt file and make it into a vector
def loadMessages(filename):
    fcorpus = open(filename, 'r')
    fcorpus = fcorpus.read() #str
    vectorCounter = makeVecs.text2vec(fcorpus)
    return (vectorCounter)



# Print cosine similarity scores
def cosineSimilarityScore(vector1, vector2):
    # cos (self, self) should = 1
    # cosine_sim_score1 = (str(makeVecs.cosine(vector1, vector1)))
    # print("cos (vec1, vec1): " + cosine_sim_score1)
    # cosine_sim_score2 = (str(makeVecs.cosine(vector2, vector2)))
    # print("cos (vec2, vec2): " + cosine_sim_score2)
    cosine_sim_score_1_2 = (makeVecs.cosine(vector1, vector2))
    # print("cos (vec1, vec2): " + cosine_sim_score_1_2)
    # print("-------------------------------------------------------")
    return cosine_sim_score_1_2


def load_corpus(corpus):
    if corpus == 'startrek':
        raw = "/home/hclent/data/corpora/startrek/105.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'darwin':
        raw = "/home/hclent/data/corpora/darwin.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'frankenstein':
        raw = "/home/hclent/data/corpora/frankenstein.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'youth':
        raw = "/home/hclent/data/corpora/youth.txt"
        corpus_vec = loadMessages(raw)
    return corpus_vec


def load_datasamples(query):
    data_samples = pickle.load(open('/home/hclent/data/18269575/data_samples_'+str(query)+'.pickle', "rb")) #pre-processed
    #print(len(data_samples)) #num docs
    data_vecs_list = loadFromDataSamples(data_samples)
    return data_vecs_list


def get_cosine_list(corpus_vec, data_vecs_list):
    cosine_list = []
    for vec_n in data_vecs_list:
        cosine_sim_score_1_2 = cosineSimilarityScore(corpus_vec, vec_n)
        score = float("{0:.4f}".format(float(cosine_sim_score_1_2)))
        cosine_list.append(score)
    return cosine_list


#take the list cosines and match scores with the url to the paper
def add_urls(query, cosine_list):
    url_list = []
    doc_list = []
    histogram_labels = [] #this is what will be in the visualization

    pmid_list = query.split('+') #list of string pmids
    for user_input in pmid_list:
        urls = db_citation_urls(user_input)
        for url in urls:
            url_list.append(url)
    num_papers = len(url_list)
    for i in range(1, (num_papers+1)):
        name = str('Doc'+str(i))
        doc_list.append(name)
    titles = list(zip(url_list, doc_list))
    for t in titles:
        label = '<a href='+str(t[0])+'>'+str(t[1])+'</a>'
        histogram_labels.append(label)
    #need to sort the cosines and return a list
    #need to sort the histogram_labels to match that order
    combo = list(zip(cosine_list, histogram_labels))
    sorted_combos = sorted(combo, reverse=False)
    return sorted_combos

def prepare_for_histogram(sorted_combos):
    x = [] #url labels
    y = [] #data points
    for combo in sorted_combos:
        x.append(combo[1]) #append url label to x
        y.append(combo[0]) #append cosine score to y
    return x, y




