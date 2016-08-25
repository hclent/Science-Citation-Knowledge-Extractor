import re
from collections import defaultdict
from multi_preprocess import *
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

py.sign_in('hclent', 'eeg49e9880')




def frequency_list(nes_list): #returns list of Dicts
    nes_stopwords = ['gene', 'genes', 'genetics', 'genome', 'chromosome', 'chromosomes', 'result', 'line', 'time']
    nesDict = defaultdict(lambda:0)
    for docs in nes_list:
        for key in docs:
            nes = (docs[key])
            for n in nes:
                if n not in nes_stopwords:
                    nesDict[n] += 1
    return nesDict


def wordcloud(nesDict, x):
    wordcloud_list = []
    for nes in nesDict:
         if int(nesDict[nes]) > x:
            #entry = {"text": nes, "size": nesDict[nes]} #no scaling
            entry = {"text": nes, "size": nesDict[nes]*.25} #scaling
            wordcloud_list.append(entry)
    wordcloud_list = re.sub('\'', "\"", str(wordcloud_list))
    return wordcloud_list



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



def plotHeatmap(x_docs, y_words, z_counts):
    data = [
        go.Heatmap(
            z = z_counts,
            x = x_docs,
            y = y_words,
        )
    ]
    #py.iplot(data, filename='labelled-heatmap')
    plot(data)
    print("plotted the heatmap!")



bdocs1 = retrieveBioDocs("18269575")
datas1, neslist1 = loadBioDoc(bdocs1)
bdocs2 = retrieveBioDocs("18952863")
datas2, neslist2 = loadBioDoc(bdocs2)

for docs in datas2:
    datas1.append(docs)
print("added datas 2 to datas 1")
for nes in neslist2:
    neslist1.append(nes)
print("added nes 2 to nes 1")



nesDict = frequency_list(neslist1)
x_docs, y_words, z_counts = doHeatmap(nesDict, 200, datas1)
plotHeatmap(x_docs, y_words, z_counts)









# biodocs2 = retrieveBioDocs("18952863")
# data_samples2, nes_list2 = loadBioDoc(biodocs2)
# for n in nes_list2:
#     nes_list1.append(n)






# Beginning and Inside
# So the *I* signals that the mention continues.
# I-Gene_or_gene_product
# B-Family
# B-CellType
# B-Simple_chemical
# B-TissueType
# I-CellType
# B-Organ
# I-TissueType
# I-Organ
# I-Cellular_component
# I-Family
# B-CellLine
# B-Gene_or_gene_product
# B-Cellular_component
# B-Species
# B-BioProcess


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

