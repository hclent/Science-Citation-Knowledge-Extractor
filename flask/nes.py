import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import re, os, json
from collections import defaultdict
from processors import *




#Make dictionary with NES counts {'gluten': 5, 'span': 9}
def frequency_dict(nes_list, category_list):
    #nes_stopwords = ['gene', 'genes', 'genetics', 'genome', 'genomic', 'chromosome', 'chromosomes', 'result', 'line', 'time']
    nesDict = defaultdict(lambda:0)
    for docs in nes_list: #docs is dict
        for key in docs: #key of dict
            for category in category_list: #no error with category 'Potato'
                if key == category:
                    nes = (docs[key])
                    for n in nes:
                        nesDict[n] += 1
                    #     if n not in nes_stopwords:
                    #        nesDict[n] += 1
    print(nesDict)
    return nesDict

#D3 wordcloud
def wordcloud(nesDict, x):
    wordcloud_list = []
    for nes in nesDict:
         if int(nesDict[nes]) > x:
            entry = {"text": nes, "size": nesDict[nes]} #no scaling
            #entry = {"text": nes, "size": nesDict[nes]*.25} #scaling
            wordcloud_list.append(entry)
    wordcloud_json = json.dumps(wordcloud_list)
    return wordcloud_json


#makes the data for plotly heatmap
def doHeatmap(nesDict, n, lemma_samples):
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

    lemma_list = [l[1] for l in lemma_samples]

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
            for document in lemma_list: #[['doc 1 words'], []]
                #docname = "doc"+str(i)
                #print(docname)

                runningDict = defaultdict(lambda:0)

                if len_word == 1:
                    for unigram in document:
                        if unigram == word:
                            #print(str(unigram)+" is the same as "+str(word))
                            runningDict[word] += 1
                        if unigram != word:
                            runningDict[word] += 0
                #bigrams
                if len_word == 2:
                    b = 2
                    bigrams = [ document[i:i+b] for i in range(len(document)-b)]
                    for grams in bigrams:
                        q = ' '
                        grams = str(q.join(grams))
                        if grams == word:
                            runningDict[word] +=1
                            #print(grams+' = '+word)
                        if grams != word:
                            runningDict[word] +=0

                #trigrams are too long for this data visualization

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
    for document in lemma_samples:
        pmcid = document[0]
        #TODO: MAKE A BETTER DOC LABEL FOR X AXIS!!!
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

#format x, y, z for Seaborn cluster map
def make_seaborn_data(x, y, z):
    seaData = []

    d = 0
    for document in x:
        i = 0
        for word in y:
            # print(i)
            count = z[i][d]  # z[i] = word, z[i][j] = count of word for document
            data = (word, document, count)
            seaData.append(data)

            i += 1
        d += 1
    return seaData

#make jpeg clusterMap
def makeClusterMap(seaData, query):
    seaFrame = pd.DataFrame((seaData), columns=['word', 'publication', 'count'])
    sea2 = seaFrame.pivot("word", "publication", "count")
    seamap = sns.clustermap(sea2)
    hm = seamap.ax_heatmap.get_position()
    plt.setp(seamap.ax_heatmap.yaxis.get_majorticklabels(), fontsize=6, rotation=0)
    seamap.ax_heatmap.set_position([hm.x0, hm.y0, hm.width, hm.height])
    col = seamap.ax_col_dendrogram.get_position()
    seamap.ax_col_dendrogram.set_position([col.x0, col.y0, col.width, col.height])

    save_path = '/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/images'  # in the folder 'clustermaps'
    save_name = ('cm_' + (str(query)) + '.png')
    completeName = os.path.join(save_path, save_name)  # with the query for a name
    seamap.savefig(completeName)
    return save_name


