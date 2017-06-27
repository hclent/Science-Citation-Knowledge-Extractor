import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import re, os, json
from collections import defaultdict
from processors import *
from database_management import db_citations_mini_hyperlink, db_citations_mini_year



#Make dictionary with NES counts {'gluten': 5, 'span': 9}
def frequency_dict(nes_list, category_list):
    nesDict = defaultdict(lambda:0)
    for docs in nes_list: #docs is dict
        for key in docs: #key of dict
            for category in category_list: #no error with category 'Potato'
                if key == category:
                    nes = (docs[key])
                    for n in nes:
                        nesDict[n] += 1
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
    years = []  # so we can sort by years
    lemma_Dict = defaultdict(lambda: 0)

    y = [] # y list of words
    x = [] # x list of documents
    z = [] # z list of lists with word counts for each document


    #STEP 1: SORT THE LEMMA_SAMPLES BY YEAR SO THAT
    #X-AXIS ON VIS CAN BE SORTED BY YEAR AS WELL :D
    for document in lemma_samples:
        pmcid = document[0]
        year = db_citations_mini_year(pmcid)
        years.append(year)


    #Sort lemma_samples by year, that way the vis goes from early years - later years
    zipped = zip(years, lemma_samples)
    sorted_data = list(sorted(zipped, key=lambda z: int(z[0])))
    sorted_lemmas = [l[1][1] for l in sorted_data] #just the lemmas in lemma_samples
    sorted_pmcids = [l[1][0] for l in sorted_data] #just the pmcids in lemma_samples

    for pmcid in sorted_pmcids:
        #print(pmcid)
        label = db_citations_mini_hyperlink(pmcid)
        keep_label = label[0] #there could be multiple, so just take the first one
        x.append(keep_label)

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
            for document in sorted_lemmas: #[['doc 1 words'], []]
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

    return x, y, z, years



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


