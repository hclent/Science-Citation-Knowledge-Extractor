import pickle, os
import naive_makeVecs as makeVecs #mine
from database_management import db_citation_urls, db_citations_hyperlink_retrieval, pmid2pmcid #mine


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


#only load eligible_papers into load_corpus
def load_corpus(corpus, eligible_papers):
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
    if corpus == 'austen':
        raw = "/home/hclent/data/corpora/austen.txt"
        corpus_vec = loadMessages(raw)
    #new
    if corpus == 'brain_speech':
        raw = "/home/hclent/data/corpora/brain_speech.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'bible':
        raw = "/home/hclent/data/corpora/bible.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'grecoroman':
        raw = "/home/hclent/data/corpora/grecoroman_med.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'last_evolution':
        raw = "/home/hclent/data/corpora/last_evolution.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'mars':
        raw = "/home/hclent/data/corpora/mars.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'mouse':
        raw = "/home/hclent/data/corpora/mouse.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'sherlock':
        raw = "/home/hclent/data/corpora/sherlock.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'yeast':
        raw = "/home/hclent/data/corpora/yeast.txt"
        corpus_vec = loadMessages(raw)
    if corpus == 'chrom1':
        raw = "/home/hclent/data/corpora/ChromNum01.txt"
        corpus_vec = loadMessages(raw)
    # For loading the query papers
    # eligible_papers = [('paper1', '18952863', '/home/hclent/data/pmcids/259/367/2593677.txt')]
    if corpus == 'paper1':
        raw = str(eligible_papers[0][2]) #'/home/hclent/data/pmcids/259/367/2593677.txt'
        corpus_vec = loadMessages(raw)
    if corpus == 'paper2':
        raw = str(eligible_papers[1][2])
        corpus_vec = loadMessages(raw)
    if corpus == 'paper3':
        raw = str(eligible_papers[2][2])
        corpus_vec = loadMessages(raw)
    if corpus == 'paper4':
        raw = str(eligible_papers[3][2])
        corpus_vec = loadMessages(raw)
    if corpus == 'paper5':
        raw = str(eligible_papers[4][2])
        corpus_vec = loadMessages(raw)
    #also gonna have to put in an exception handler for if there is no PMCID for the input paper (thus no corpus)
    #can return a flash error for that
    return corpus_vec


def load_datasamples(query):
    data_samples = pickle.load(open('/home/hclent/data/data_samples/data_samples_'+str(query)+'.pickle', "rb")) #pre-processed
    #print(len(data_samples)) #num docs
    data_vecs_list = loadFromDataSamples(data_samples)
    return data_vecs_list


def get_cosine_list(corpus_vec, data_vecs_list):
    cosine_list = []
    for vec_n in data_vecs_list:
        cosine_sim_score_1_2 = cosineSimilarityScore(corpus_vec, vec_n)
        score = float("{0:.4f}".format(float(cosine_sim_score_1_2)))
        score = float(score * 100) #.25 --> 25%
        cosine_list.append(score)
    return cosine_list


#take the list cosines and match scores with the url to the paper
def add_urls(query, cosine_list):
    url_list = []
    doc_list = []
    histogram_labels = [] #this is what will be in the visualization
    apa_labels = []

    pmid_list = query.split('+') #list of string pmids
    for user_input in pmid_list:
        #get the urls
        urls = db_citation_urls(user_input)
        for url in urls:
            url_list.append(url)
        #get hyperlinked apa citations for click event
        citations = db_citations_hyperlink_retrieval(user_input)
        for c in citations:
            apa_labels.append(c)

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
    combo = list(zip(cosine_list, histogram_labels, apa_labels))
    sorted_combos = sorted(combo, reverse=False)
    return sorted_combos

def prepare_for_histogram(sorted_combos):
    x = [] #url labels
    y = [] #data points
    names = []
    for combo in sorted_combos:
        x.append(combo[1]) #append url label to x
        y.append(combo[0]) #append cosine score to y
        names.append(combo[2])
    return x, y, names


