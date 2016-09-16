import makeVecs as makeVecs
import math, string, re
from collections import Counter


#Input: publication
#Output: normalized vector (counter), representing stringIn.

#Example
doc1 = "/home/hclent/data/18269575/18269575_62.txt"
doc2 = "/home/hclent/data/18269575/18269575_63.txt"
scifi = "/home/hclent/data/corpora/startrek/107.txt"


def loadMessages(filename):
    fcorpus = open(filename, 'r')
    fcorpus = fcorpus.read() #str

    vectorCounter = makeVecs.text2vec(fcorpus)

    return (vectorCounter)


vecs1 = loadMessages(doc1)
print(vecs1)
vecs2 = loadMessages(doc2)
print("################################################################################################")
vecs3= loadMessages(scifi)
print(vecs3)


#Print cosine similarity scores
def cosineSimilarityScore(vector1, vector2, vector3):

    cosine_sim_score1 = (str(makeVecs.cosine(vector1, vector1)))
    print("cos (vec1, vec1): " + cosine_sim_score1)


    cosine_sim_score2 = (str(makeVecs.cosine(vector2, vector2)))
    print("cos (vec2, vec2): " + cosine_sim_score2)


    cosine_sim_score_1_2 = (str(makeVecs.cosine(vector1, vector2)))
    print("cos (vec1, vec2): " + cosine_sim_score_1_2)


    cosine_sim_score_1_3 = (str(makeVecs.cosine(vector1, vector3)))
    print("cos (vec1, vec3): " + cosine_sim_score_1_3)


    cosine_sim_score_2_3 = (str(makeVecs.cosine(vector2, vector3)))
    print("cos (vec2, vec3): " + cosine_sim_score_2_3)


cosineSimilarityScore(vecs1, vecs2, vecs3)



############### EXAMPLE #################
# cos (vec1, vec1): 0.9999999999999983 #self
# cos (vec2, vec2): 1.0000000000000064 #self
# cos (vec1, vec2): 0.21345131811302565 #two scientific texts
# cos (vec1, vec3): 0.03610515107865741 #scientific versus movie review
# cos (vec2, vec3): 0.03499071664127394 #scientific versus movie review