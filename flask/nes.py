import re
from collections import defaultdict
from processors import *
from multi_preprocess import * #mine #connects to NLP Server


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

                #trigrams
                if len_word == 3:
                    t = 3
                    trigrams = [unigrams_list[i:i+t] for i in range(len(unigrams_list)-t)]
                    for grams in trigrams:
                        q = ' '
                        grams = str(q.join(grams))
                        if grams == word:
                            runningDict[word] +=1
                            #print(grams+' = '+word)
                        if grams != word:
                            runningDict[word] +=0

                #4 grams
                if len_word > 3:
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


# def get_data_and_ner(pmid):
# 	biodocs = retrieveBioDocs(str(pmid)) #a bunch of strings
# 	data_samples, neslist = loadBioDoc(biodocs)
# 	return data_samples, neslist


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
# category_list = ['BioProcess', 'CellLine', 'Cellular_component', 'Family', 'Gene_or_gene_product', 'Organ', 'Simple_chemical', 'Site', 'Species', 'TissueType']
# nesDict = frequency_dict(neslist, category_list)
# x, y, z = doHeatmap(nesDict, 100, data_samples)
#
# print("X:")
# print(x)
# print("######################################################################")
# print("Y:")
# print(y)
# print("######################################################################")
# print("Z:")
# print(z)
# print("######################################################################")