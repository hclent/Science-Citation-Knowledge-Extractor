import naive_makeVecs as makeVecs
import pickle
from database_management import db_citation_titles


#Load from pickled data_samples instead of filename
def loadFromDataSamples(data_samples):
    vecs_list = []

    for document in data_samples:
        vectorCounter = makeVecs.text2vec(document)
        #print(vectorCounter)
        vecs_list.append(vectorCounter)
    return vecs_list


#Load txt file
def loadMessages(filename):
    fcorpus = open(filename, 'r')
    fcorpus = fcorpus.read() #str

    vectorCounter = makeVecs.text2vec(fcorpus)

    return (vectorCounter)



# Print cosine similarity scores
def cosineSimilarityScore(vector1, vector2):

    cosine_sim_score1 = (str(makeVecs.cosine(vector1, vector1)))
    print("cos (vec1, vec1): " + cosine_sim_score1)


    cosine_sim_score2 = (str(makeVecs.cosine(vector2, vector2)))
    print("cos (vec2, vec2): " + cosine_sim_score2)


    cosine_sim_score_1_2 = (str(makeVecs.cosine(vector1, vector2)))
    print("cos (vec1, vec2): " + cosine_sim_score_1_2)
    print("-------------------------------------------------------")
    return cosine_sim_score_1_2



star_trek = "/home/hclent/data/corpora/startrek/105.txt"
vecs1 = loadMessages(star_trek)


data_samples = pickle.load(open("/home/hclent/data/18269575/data_samples_18952863+18269575.pickle", "rb")) #pre-processed
print(len(data_samples))

vecs_list = loadFromDataSamples(data_samples)
cosine_list = []
for vec_n in vecs_list:
    cosine_sim_score_1_2 = cosineSimilarityScore(vecs1, vec_n)
    cosine_list.append(cosine_sim_score_1_2)


print(cosine_list)

doc_list = []

for i in range(1,166):
    name = str('Doc'+str(i))
    doc_list.append(name)
    i+=1

print(doc_list)


hover = []

titles= db_citation_titles(18269575)
for t in titles:
    hover.append(t)
    
titles = db_citation_titles(18952863)
for t in titles:
    hover.append(t)




print(hover)
