import re
from collections import defaultdict
from processors import *
import numpy as np
from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans
import pickle, random, time
from multi_preprocess import retrieveBioDocs, loadBioDoc



#Make dictionary with NES counts {'gluten': 5, 'span': 9}
def frequency_dict(nes_list, category_list):
    #nes_stopwords = ['gene', 'genes', 'genetics', 'genome', 'genomic', 'chromosome', 'chromosomes', 'result', 'line', 'time']
    nesDict = defaultdict(lambda:0)
    for docs in nes_list:
        for key in docs:
            for category in category_list: #no error with category 'Potato'
                if key == category:
                    nes = (docs[key])
                    for n in nes:
                        nesDict[n] += 1
                    #     if n not in nes_stopwords:
                    #        nesDict[n] += 1
    return nesDict

#D3 wordcloud
def wordcloud(nesDict, x):
    wordcloud_list = []
    for nes in nesDict:
         if int(nesDict[nes]) > x:
            entry = {"text": nes, "size": nesDict[nes]} #no scaling
            #entry = {"text": nes, "size": nesDict[nes]*.25} #scaling
            wordcloud_list.append(entry)
    wordcloud_list = re.sub('\'', "\"", str(wordcloud_list))
    return wordcloud_list


#makes the data for plotly heatmap
def doHeatmap(nesDict, n, data_samples):
    y = [] # y list of words
    x = [] # x list of documents
    z = [] # z list of lists with word counts for each document

    lemma_Dict = defaultdict(lambda:0)


    #remove plurals (very rough)
    for word in nesDict:
        p1 = re.compile('[a-z]{2,}[b-df-hj-np-tv-z]{1,}(es$|s$)')
        match1 = p1.search(word)
        if not match1:
            lemma_Dict[word] += nesDict[word]


    for word in lemma_Dict:
        # print(word)
        split_word = word.split(" ")
        # print(split_word)
        # print(type(split_word))
        len_word = len(split_word)
        # print(len_word)
        if int(lemma_Dict[word]) > int(n):
            y.append(word)
            #i = 0
            #print(word)
            word_counts = []
            for documents in data_samples:
                #docname = "doc"+str(i)
                #print(docname)
                unigrams_list = documents.split(" ")
                runningDict = defaultdict(lambda:0)

                if len_word == 1:
                    for unigram in unigrams_list:
                        if unigram == word:
                            #print(str(unigram)+" is the same as "+str(word))
                            runningDict[word] += 1
                        if unigram != word:
                            runningDict[word] += 0
                #bigrams
                if len_word == 2:
                    b = 2
                    bigrams = [ unigrams_list[i:i+b] for i in range(len(unigrams_list)-b)]
                    for grams in bigrams:
                        q = ' '
                        grams = str(q.join(grams))
                        if grams == word:
                            runningDict[word] +=1
                            #print(grams+' = '+word)
                        if grams != word:
                            runningDict[word] +=0

                #trigrams are too long for this data visualization
                # if len_word == 3:
                #     t = 3
                #     trigrams = [unigrams_list[i:i+t] for i in range(len(unigrams_list)-t)]
                #     for grams in trigrams:
                #         q = ' '
                #         grams = str(q.join(grams))
                #         if grams == word:
                #             runningDict[word] +=1
                #             #print(grams+' = '+word)
                #         if grams != word:
                #             runningDict[word] +=0

                #3 grams, 4 grams, +
                if len_word > 2:
                    #print("the named entity is too long")
                    pass
                #print(runningDict)
                word_counts.append(runningDict[word])
                #i +=1
            #print(word_counts)
            z.append(word_counts)

    i = 1
    for documents in data_samples:
        docname = "doc"+str(i)
        x.append(docname)
        i += 1
    return x, y, z


#deciding k for num of clusters
# def elbowMethod(z):
#     z_array = np.array(z)
#     linkArray =  linkage(z_array, 'average')
#     last = linkArray[-10:, 2]
#     # last_rev = last[::-1]
#     # idxs = np.arange(1, len(last) + 1)
#     acceleration = np.diff(last, 2)  # 2nd derivative of the distances
#     acceleration_rev = acceleration[::-1]
#     k = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters
#     if k > 6:
#         k = 6
#         print("k is greater than 6 yikes")
#     else:
#         print("clusters: " + str(k))
#     return k


#perform k-means on y-axis values (i.e. z matrix)
def do_kmeans_on_y(z, k_clusters):
    sparse_matrix = np.array(z)
    t0 = time.time()
    #print("* Beginning k-means clustering ... ")
    num_clusters = int(k_clusters)
    km = KMeans(init='k-means++', n_clusters=num_clusters)
    km.fit(sparse_matrix)
    clusters = km.labels_.tolist()
    #print("done in %0.3fs." % (time.time() - t0))
    return clusters

#re-order the y-axis now for the new clustered heatmap
def sort_y(clusters, y):
    pairs = list(zip(clusters, y))
    # print(pairs)
    list0 = []
    list1 = []
    list2 = []
    list3 = []
    list4 = []
    list5 = []

    for p in pairs:
        if p[0] == 0:
            list0.append(p[1])
        if p[0] == 1:
            list1.append(p[1])
        if p[0] == 2:
            list2.append(p[1])
        if p[0] == 3:
            list3.append(p[1])
        if p[0] == 4:
            list4.append(p[1])
        if p[0] == 5:
            list4.append(p[1])
    new_order = list0 + list1 + list2 + list3 + list4 + list5
    return new_order


# heatmap coordinates for handling a clustered y-axis
def doClustermap(nesDict, n, data_samples, new_order, old_x):
    x = []
    y = []
    z = []

    lemma_Dict = defaultdict(lambda: 0)
    for word in nesDict:
        p1 = re.compile('[a-z]{2,}[b-df-hj-np-tv-z]{1,}(es$|s$)')
        match1 = p1.search(word)
        if not match1:
            lemma_Dict[word] += nesDict[word]
    #print("* made lemmas dict ...")

    for word in new_order:

        # for word in lemma_Dict:
        split_word = word.split(" ")
        len_word = len(split_word)
        if int(lemma_Dict[word]) > int(n):
            y.append(word)
            word_counts = []
            for j in range(0, len(data_samples)):  # go in order. nes are clustered, not docs
                documents = data_samples[int(j)]
                unigrams_list = documents.split(" ")
                runningDict = defaultdict(lambda: 0)
                if len_word == 1:
                    for unigram in unigrams_list:
                        if unigram == word:
                            runningDict[word] += 1
                        if unigram != word:
                            runningDict[word] += 0
                if len_word == 2:
                    b = 2
                    bigrams = [unigrams_list[i:i + b] for i in range(len(unigrams_list) - b)]
                    for grams in bigrams:
                        q = ' '
                        grams = str(q.join(grams))
                        if grams == word:
                            runningDict[word] += 1
                            # print(grams+' = '+word)
                        if grams != word:
                            runningDict[word] += 0
                if len_word > 2:
                    pass
                word_counts.append(runningDict[word])
            z.append(word_counts)

    for i in range(1, len(old_x) + 1):
        docname = "doc" + str(i)
        x.append(docname)
    #print("* done with stuff ...")
    return x, y, z




def get_data_and_ner(pmid):
    biodocs = retrieveBioDocs(str(pmid)) #a bunch of strings
    #print(biodocs)
    data_samples, neslist, sents, tokens = loadBioDoc(biodocs)
    return data_samples, neslist


# data_samples = []
# neslist = []
#
# ds1, nes1 = get_data_and_ner(18269575)
# ds2, nes2 = get_data_and_ner(18952863)
#
# for d in ds1:
#     data_samples.append(d)
# for d in ds2:
#     data_samples.append(d)
# for n in nes1:
#     neslist.append(n)
# for n in nes2:
#     neslist.append(n)
# print("put all the things in the lists lol")
#
#
# #category_list = ['BioProcess', 'CellLine', 'Cellular_component', 'Family', 'Gene_or_gene_product', 'Organ', 'Simple_chemical', 'Site', 'Species', 'TissueType']
# category_list = ['Species']
#
# nesDict = frequency_dict(neslist, category_list)
# x, y, z = doHeatmap(nesDict, 10, data_samples)
#
# k =  elbowMethod(z)
# clusters = do_kmeans_on_y(z, 6)
# new_order = sort_y(clusters, y)
# x1, y1, z1 = doClustermap(nesDict, 10, data_samples, new_order, x)
# print(y1)





# print("X:")
# print(x)
# print("######################################################################")
# print("Y:")
# print(y)
# print("######################################################################")
# print("Z:")
# print(z)
# print("######################################################################")