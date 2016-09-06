import re
from collections import defaultdict
#import plotly.plotly as py
#import plotly.graph_objs as go
#from plotly.offline import plot
#py.sign_in('hclent', 'eeg49e9880')
from multi_preprocess import * #mine #connects to NLP Server


###### changes to py-processors NLP #########


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

    for word in nesDict:
        if nesDict[word] > n:
            y.append(word)
            #i = 0
            #print(word)
            word_counts = []
            for documents in data_samples:
                #docname = "doc"+str(i)
                #print(docname)
                unigrams_list = documents.split(" ")
                runningDict = defaultdict(lambda:0)
                for unigram in unigrams_list:
                    if unigram == word:
                        #print(str(unigram)+" is the same as "+str(word))
                        runningDict[word] += 1
                    if unigram != word:
                        runningDict[word] += 0
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

#uses doHeatmap data and generates plotly graph
def plotHeatmap(x_docs, y_words, z_counts):
    data = [
        go.Heatmap(
            z = z_counts,
            x = x_docs,
            y = y_words,
            colorscale=[[0.0, 'rgb(204,204,204)'], [0.01111111111111111, 'rgb(69,79,220)'],[0.1111111111111111, 'rgb(84,73,210)'], [0.2222222222222222, 'rgb(98,67,201)'], [0.3333333333333333, 'rgb(127,54,181)'], [0.4444444444444444, 'rgb(142,48,172)'], [0.5555555555555556, 'rgb(171,36,152)'], [0.6666666666666666, 'rgb(199,24,133)'], [0.7777777777777778, 'rgb(214,17,123)'], [0.8888888888888888, 'rgb(228,11,114)'], [1.0, 'rgb(243,5,104)']]
        )
    ]
    #py.iplot(data, filename='labelled-heatmap') #plot to online
    plot(data,  filename='blahblahblah.html')
    print("plotted the heatmap!")



# bdocs1 = retrieveBioDocs("18269575") #a bunch of strings
# data_samples, neslist1 = loadBioDoc(bdocs1)
# chosen_categories = ['CellLine', 'Cellular_component','Gene_or_gene_product', 'Site']
# frequency_dict(neslist1, chosen_categories)

# wcl = wordcloud(nesDict, 50)
# print(wcl)


# x_docs, y_words, z_counts = doHeatmap(nesDict, 100, data_samples)
# plotHeatmap(x_docs, y_words, z_counts)
#sortCategories(data_samples, neslist1)



########## Graveyard ###########
#d3
#{"data": [ {"timestamp": "time", "value": {"PM2.5": 30}}, ... , ...   ]     }
#{"data": [{"document": "doc1", "value": {"tissue": 30}}, ... , ...] }
# def heatmap_vis(data_samples, frequency_list):
#     hmDict = {"data": []}
#
#     i = 0
#     for document in data_samples:
#         docname = "doc"+str(i)
#         #runningDict = {"document": docname, "value": defaultdict(lambda:0)}
#         text = document.split(' ') #split on white space
#         for nes in frequency_list: #nes is a dictionary entry
#             word = str(nes["text"])
#             x = defaultdict(lambda:0)
#             for t in text:
#                 if word == t:
#                     #print(word + " is in "+ str(docname) )
#                     x[word] +=1
#                 else:
#                     x[word] += 0
#             runningDict = {"document": docname, "value": dict(x) }
#             #print(runningDict)
#             (hmDict["data"]).append(runningDict)
#         i += 1
#     print(hmDict)

