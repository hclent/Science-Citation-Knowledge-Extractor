import math
from collections import Counter
import string
import nltk
from nltk.corpus import stopwords


#Cosine similarity
#Input: Two vectors
#Output: Dot product of the vectors
def cosine(vec1, vec2):
    return dotProduct(vec1, vec2)


#Get the dot product
#Input two vectors (type Counters)
def dotProduct(vec1, vec2):
    runningSum = float(0.0)
    for val in vec1:
        value1 = vec1[val]
        value2 = vec2[val]
        prod = float(value1) * float(value2)
        runningSum += prod

    return runningSum


#Calculate vector length for normalization
def vecLength(vec):
    runningSum = float(0.0)
    for value in vec:
        exp = vec[value] ** 2
        runningSum += exp

    length = math.sqrt(runningSum)
    return length


#normalize
#Input: Vector
#Output: Normalized vector
def normalize(vec):
    length = vecLength(vec)
    for val in vec:
        vec[val] = vec[val] / length


#Convert text to a vector space
#Input: Text
#Output: Vector
def text2vec(stringIn):
    #remove punctuation
    translator = str.maketrans({key: None for key in string.punctuation}) #Updated translator for python 3.x
    filteredString1 = (stringIn.translate(translator))

    #convert to lower-case
    filteredString2 = filteredString1.lower()

    #split to words
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