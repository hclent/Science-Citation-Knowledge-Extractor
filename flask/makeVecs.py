import math
from collections import Counter
import string
import nltk
from nltk.corpus import stopwords

#Vector functions for cosine similarity based classification
#Vectors must be normalized in order for this to work!


#cosine similarity
def cosine(vec1, vec2):
    return dotProduct(vec1, vec2)


#get the dot product
def dotProduct(vec1, vec2):
    runningSum = float(0.0)
    for key in vec1:
        value1 = vec1[key]
        value2 = vec2[key]
        prod = float(value1) * float(value2)
        runningSum += prod

    return runningSum


#calculate vector length for normalization
def vecLength(vec):
    runningSum = float(0.0)
    for key in vec:
        exp = vec[key] ** 2
        runningSum += exp

    length = math.sqrt(runningSum)
    return length


#normalize
def normalize(vec):
    length = vecLength(vec)
    for key in vec:
        vec[key] = vec[key] / length



def text2vec(stringIn):
    #remove punctuation
    translator = str.maketrans({key: None for key in string.punctuation}) #Updated translator for python 3.x
    filteredString1 = (stringIn.translate(translator))

    #convert to lower-case
    filteredString2 = filteredString1.lower()

    #pplit to words
    words = filteredString2.split()

    #must delete stop words for cosine similarity scores to be representative of document's content!
    stopwords = nltk.corpus.stopwords.words('english')
    words = [w for w in words if w.lower() not in stopwords]

    #convert to vector
    outVec = Counter()
    for word in words:
        outVec[word] += 1


    #normalize
    normalize(outVec)

    return outVec
