import string, os
import naive_makeVecs as makeVecs #mine
from database_management import db_citations_mini_hyperlink, db_citations_hyperlink_retrieval, db_pmid_axis_label, db_pmid_hyperlink_retrieval #mine
from cache_lemma_nes import load_lemma_cache #mine
from configapp import app


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
    cosine_sim_score_1_2 = (makeVecs.cosine(vector1, vector2))
    return cosine_sim_score_1_2


#only load eligible_papers into load_corpus
def load_corpus(corpus, eligible_papers):
    if corpus == 'startrek':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), 'startrek.txt')
        corpus_vec = loadMessages(raw)
        color = 'rgb(63, 100, 168)'
    if corpus == 'darwin':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), 'darwin.txt')
        corpus_vec = loadMessages(raw)
        color = 'rgb(8, 114, 32)'
    if corpus == 'frankenstein':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "frankenstein.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(92, 59, 107)'
    if corpus == 'youth':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "youth.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(63, 100, 168)'
    if corpus == 'austen':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "austen.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(191, 110, 167)'
    #new
    if corpus == 'brain_speech':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "brain_speech.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(8, 114, 32)'
    if corpus == 'bible':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "bible.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(92, 59, 107)'
    if corpus == 'grecoroman':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "grecoroman_med.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(8, 114, 32)'
    if corpus == 'last_evolution':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "last_evolution.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(63, 100, 168)'
    if corpus == 'mars':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "mars.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(63, 100, 168)'
    if corpus == 'mouse':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "mouse.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(8, 114, 32)'
    if corpus == 'sherlock':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "sherlock.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(92, 59, 107)'
    if corpus == 'yeast':
        raw = os.path.join((app.config['PATH_TO_CORPORA']), "yeast.txt")
        corpus_vec = loadMessages(raw)
        color = 'rgb(8, 114, 32)'
    # For loading the query papers
    # eligible_papers = [('paper1', '18952863', '/home/hclent/data/pmcids/259/367/2593677.txt')]
    if corpus == 'paper1':
        raw = str(eligible_papers[0][2]) #'~/data/pmcids/259/367/2593677.txt'
        corpus_vec = loadMessages(raw)
        color = 'rgb(8, 114, 32)'
    if corpus == 'paper2':
        raw = str(eligible_papers[1][2])
        corpus_vec = loadMessages(raw)
        color = 'rgb(8, 114, 32)'
    if corpus == 'paper3':
        raw = str(eligible_papers[2][2])
        corpus_vec = loadMessages(raw)
        color = 'rgb(8, 114, 32)'
    if corpus == 'paper4':
        raw = str(eligible_papers[3][2])
        corpus_vec = loadMessages(raw)
        color = 'rgb(8, 114, 32)'
    if corpus == 'paper5':
        raw = str(eligible_papers[4][2])
        corpus_vec = loadMessages(raw)
        color = 'rgb(8, 114, 32)'
    #also gonna have to put in an exception handler for if there is no PMCID for the input paper (thus no corpus)
    #can return a flash error for that
    return corpus_vec, color


#Updated to use new cache :) yay
def load_datasamples(query):
    lemma_samples = load_lemma_cache(query)
    lemma_list = [l[1] for l in lemma_samples]
    pmcids_list = [l[0] for l in lemma_samples]
    list_data_strings = [' '.join(map(str, l)) for l in  lemma_list]
    data_vecs_list = loadFromDataSamples(list_data_strings)
    return data_vecs_list, pmcids_list


def get_cosine_list(corpus_vec, data_vecs_list):
    cosine_list = []
    for vec_n in data_vecs_list:
        cosine_sim_score_1_2 = cosineSimilarityScore(corpus_vec, vec_n)
        score = float("{0:.4f}".format(float(cosine_sim_score_1_2)))
        score = float(score * 100) #.25 --> 25%
        cosine_list.append(score)
    return cosine_list


#get cosine scores for any eligible papers
def get_cosine_eligible(corpus_vec, eligible_papers):
    eligible_cosines = []
    if len(eligible_papers) > 0:
        for paper in eligible_papers:
            filename = paper[2]
            vector_counter = loadMessages(filename)
            cosine_score = cosineSimilarityScore(corpus_vec, vector_counter)
            score = float("{0:.4f}".format(float(cosine_score)))
            score = float(score * 100)  # .25 --> 25%
            eligible_cosines.append(score)
    return eligible_cosines





#take the list cosines and match scores with the url to the paper
#Updated with new titles, making sure there's no repeats :)
def add_urls(cosine_list, color, pmcids_list, conn):

    histogram_labels = [] #this is what will be in the visualization
    apa_labels = []

    alphabet = list(string.ascii_lowercase)
    for pmcid in pmcids_list:
        label = db_citations_mini_hyperlink(pmcid, conn)
        keep_label = label[0]  # there could be multiple records from db, so just take the first one
        # Step 1: check if its in the x list
        repeat_count = histogram_labels.count(keep_label)
        if repeat_count > 0:
            # eww hacky yucky i'm really sorry!
            add_letter = alphabet[repeat_count]
            keep_label = keep_label[:4] + add_letter + keep_label[4:]
        elif repeat_count == 0:
            pass
        histogram_labels.append(keep_label)
        #get the hyperlink apa_lables
        hyperlink_list = db_citations_hyperlink_retrieval(pmcid, conn)
        keep_hyperlink = hyperlink_list[0]
        apa_labels.append(keep_hyperlink)

    colors_list = [color] * int(len(histogram_labels))

    #need to sort the histogram_labels to match that order
    combo = list(zip(cosine_list, histogram_labels, apa_labels, colors_list))
    sorted_combos = sorted(combo, reverse=False)
    return sorted_combos



#eligible_papers = [('paper1', '18952863', '/home/hclent/data/pmcids/259/367/2593677.txt')]
def add_eligible_cosines(sorted_combos, eligible_papers, eligible_cosines, conn):
    labels_thus_far = [c[1] for c in sorted_combos]

    alphabet = list(string.ascii_lowercase)
    #get the label for the eligible papers!

    if len(eligible_papers) > 0:
        histogram_labels = []
        click_label = []
        for e in eligible_papers:
            pmid = e[1]
            label = db_pmid_axis_label(pmid, conn)
            keep_label = label[0] #there might be multiple, just get first

            #check if the label is unique from others. plotly deletes duplicate lables :/
            repeat_count = labels_thus_far.count(keep_label)
            if repeat_count > 0:
                add_letter = alphabet[repeat_count]
                keep_label = keep_label[:4] + add_letter + keep_label[4:]
            elif repeat_count == 0:
                pass

            hyperlink = db_pmid_hyperlink_retrieval(pmid, conn)
            keep_hyperlink = hyperlink[0]

            histogram_labels.append(keep_label)

            click_label.append(keep_hyperlink)

        color = 'rgb(244, 241, 48)' #Hilight inputPapers in Yellow :)
        colors_list = [color] * int(len(histogram_labels))
        combo = list(zip(eligible_cosines, histogram_labels, click_label, colors_list))
        all_combos = sorted_combos + combo
        sorted_all_combos =  sorted(all_combos, reverse=False)
        return sorted_all_combos

    #if there are no eligible papers just pass
    if len(eligible_papers) == 0:
        return sorted_combos


def prepare_for_histogram(sorted_combos):
    x = [] #url labels
    y = [] #data points
    names = []
    color = []
    for combo in sorted_combos:
        x.append(combo[1]) #append url label to x
        y.append(combo[0]) #append cosine score to y
        names.append(combo[2])
        color.append(combo[3])
    return x, y, names, color
