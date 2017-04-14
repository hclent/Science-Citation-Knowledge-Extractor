import re, os.path, sys
from scipy import spatial
import numpy as np
from gensim.models import Word2Vec
from fasttext import load_model

#make similarity matrix to decide what is/isn't connected in force dir graph
#words in same topic are connected (1.0)
#words in different topics are connected if cos(w1,w2) > some x
#words are assigned a group id for the topic they fall in
#Sorry this is really ugly :'(
def make_matrix(results, model):
    val_matrix = []

    for tup1 in results:
        topic1 = tup1[0]
        word1 = tup1[1][0]
        word1count = tup1[1][1]

        temp_row = []
        for tup2 in results:
            topic2 = tup2[0]
            word2 = tup2[1][0]

            if topic1 == topic2: #if topic number = topic number
               # print(str(word1) + " == " + str(word2))
                if word1 == word2: # if word = word
                    temp_row.append(word1count) #append word count
                if word1 != word2: #if they are not the same word, but in the same group
                    temp_row.append(1.0) #just append 1.0

            if topic1 != topic2: #if the topic numbers are different
                #print("cos(" + str(word1) + ", " + str(word2) + ")")
                #PROBLEM: Need the average embedding for the noun phrase.
                words_one = word1.split(' ') #list of words in word1
                num_words_one = len(words_one) #number of words in word1
                words_two = word2.split(' ') #list of words in word2
                num_words_two = len(words_two) #number of words in word2

                zeroes = [0.0] * 100
                sum1 = np.asarray(zeroes).T  # init empty column vector (100,)
                sum2 = np.asarray(zeroes).T  # init empty column vector (100,)
                for w1s in words_one:
                    try:
                        np_vec = model[w1s]
                        sum1 += np.add(sum1, np_vec)
                    except Exception as e1:
                        pass
                for w2s in words_two:
                    try:
                        np_vec = model[w2s]
                        sum2 += np.add(sum2, np_vec)
                    except Exception as e2:
                        pass
                #get the average embedding for the noun phrase
                avg_w1_vec = np.divide(sum1, num_words_one).T.tolist() #make a row vector, then make it a list
                avg_w2_vec = np.divide(sum2, num_words_two).T.tolist()
                #take the cosine
                cos_sim = 1 - spatial.distance.cosine(avg_w1_vec, avg_w2_vec)
                if cos_sim == 'nan':
                    #print("NANANANANNANA :( ")
                    pass
                temp_row.append(cos_sim)

        #group val is last val on each row
        group_num = int(topic1)
        temp_row.append(group_num)

        #print(temp_row)
        val_matrix.append(temp_row)

    #add final rows column
    last_row_vals = [ int(tup[0] ) for tup in results ] #this is the topic ids
    #print(last_row_vals)
    val_matrix.append(last_row_vals)

    return val_matrix


def make_csv(val_matrix, results, query):
    words_list = [words[1][0] for words in results]
    formatted_words_list = [ re.sub("\s", "_", w) for w in words_list]

    top_row = str(str(",")+(",".join(formatted_words_list)) + str(",group"))

    save_path = '/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/csvgraphs/'  # in the folder of the last pmid
    completeName = os.path.join(save_path, ('fgraph_' + (str(query)) + '.csv'))  # with the query for a name
    sys.stdout = open(completeName, "w")

    print(top_row)
    i = 0
    rows = formatted_words_list + ["group"]
    for label in rows: #label is str
        values = val_matrix[i] #list
        line = str( label + "," + str(values))
        line = re.sub( "\s+" , "", line)
        line = re.sub("\]|\[", "", line)
        line = re.sub(',nan', ',0.00001', line) #with scipy some items are nans
        print(line)
        i+=1


# results = [(2, ('gene family', 525)), (5, ('a. thaliana', 409)), (6, ('genome duplication', 282)), (2, ('gene duplication', 267)), (2, ('gene pair', 245)), (3, ('arabidopsis thaliana', 227)), (8, ('plant species', 224)), (2, ('gene loss', 203)), (9, ('expression pattern', 203)), (6, ('duplication event', 194)), (6, ('genome sequence', 180)), (1, ('protein sequence', 168)), (3, ('synteny block', 163)), (8, ('plant genome', 153)), (9, ('expression level', 152)), (0, ('wgd event', 147)), (2, ('gene expression', 143)), (5, ('b. rapa', 140)), (9, ('expression profile', 138)), (0, ('linkage group', 135)), (1, ('sequence similarity', 135)), (0, ('bac clone', 131)), (7, ('amino acid', 124)), (2, ('gene order', 119)), (6, ('reference genome', 116)), (4, ('m. truncatulum', 116)), (2, ('gene model', 111)), (0, ('mhc class', 107)), (6, ('genome assembly', 107)), (4, ('g. max', 106)), (4, ('\\ u2003', 106)), (3, ('synteny analysis', 100)), (9, ('transcription factor', 99)), (0, ('datum set', 99)), (0, ('mam gene', 98)), (1, ('sequence alignment', 98)), (2, ('gene cluster', 98)), (4, ('t \\ t', 96)), (2, ('protein coding gene', 92)), (5, ('c. arabica', 92)), (0, ('table s1', 89)), (8, ('land plant', 82)), (4, ('k value', 82)), (2, ('gene annotation', 80)), (2, ('gene copy', 79)), (2, ('target gene', 79)), (0, ('arabidopsis lyra', 79)), (8, ('flowering plant', 78)), (5, ('b cell', 78)), (4, ('l. angustifolius', 77)), (0, ('fad3 gene', 76)), (0, ('qtl region', 69)), (0, ('pcr product', 69)), (6, ('genome size', 66)), (2, ('candidate gene', 66)), (0, ('divergence time', 65)), (0, ('prr gene', 65)), (5, ('c. canephora', 65)), (4, ('a. lyra', 64)), (4, ('l. japonicus', 64)), (9, ('stress response', 62)), (0, ('dna sequence', 61)), (0, ('brassica rapa', 61)), (4, ('p. vulgari', 61)), (0, ('default parameter', 61)), (2, ('gene content', 61)), (8, ('flowering time', 60)), (3, ('maize genome', 59)), (0, ('iaa gene', 58)), (0, ('branch length', 58)), (0, ('dart marker', 58)), (8, ('model plant', 58)), (7, ('amino acid sequence', 57)), (3, ('arabidopsis genome', 57)), (0, ('ssr marker', 57)), (9, ('transcript abundance', 55)), (6, ('genome annotation', 55)), (6, ('genome evolution', 55)), (2, ('copy number', 55)), (2, ('gene tree', 54)), (2, ('gene function', 54)), (3, ('grape genome', 54)), (3, ('blast search', 53)), (9, ('transcript level', 53)), (4, ('m. itineran', 53)), (0, ('core eudicot', 53)), (0, ('hormone treatment', 52)), (0, ('acid qtl', 52)), (2, ('gene sequence', 52)), (8, ('plant lineage', 52)), (5, ('b. distachyon', 51)), (8, ('legume species', 51)), (3, ('grape ap gene', 51)), (0, ('ap gene', 51)), (0, ('table s2', 51)), (6, ('tef genome', 50)), (2, ('gene set', 50)), (0, ('reading frame', 50)), (6, ('genome duplication event', 50)), (2, ('gene structure', 50))]
# model = load_model('/home/hclent/tmp/fastText/16kmodel.vec')
# print("model loaded!")
# val_matrix = make_matrix(results, model)
# make_csv(val_matrix, results)