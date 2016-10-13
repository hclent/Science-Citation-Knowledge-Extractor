import gensim, logging, os, codecs, math
from gensim.models import Word2Vec
import numpy as np
from numpy import linalg as LA


logging.basicConfig(filename='.2vec.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
logging.info('Started')


#lacks pre-processing :(
class MySentences(object):
    def __init__(self, dirname):
        self.dirname = dirname

    def __iter__(self):
        for fname in os.listdir(self.dirname):
            if ".txt" in fname:
                for line in codecs.open(os.path.join(self.dirname, fname), "r", encoding='utf-8', errors='ignore'):
                    if len(line.split()) != 0: #make sure its not empty
                        #yield line.split()  #or line.lower().split() for lowercase
                        yield line.lower().split()


def create_model(path_to_sentences, model_name):
    print("* retrieving sentences ... ")
    sentences = MySentences(path_to_sentences) # a memory-friendly iterator
    print("* successfully retrieved sentences !!!")
    print("* making the model ... ")
    model = gensim.models.Word2Vec(sentences, min_count=1)
    print("* successfully made the model !!!")
    print("* saving the model ... ")

    #more_sentences = MySentences('/home/hclent/data/18952863')
    #model.train(more_sentences)

    model.save(str(path_to_sentences+model_name))
    print("* successfully saved the model !!!")

#create_model('/home/hclent/data/corpora/startrek/', 'startrek_model')
#create_model('/home/hclent/data/18269575', 'coge_model')

def load_model(path_to_sentences, model_name):
    print("* loading the model ... ")
    model = Word2Vec.load(str(path_to_sentences+model_name))
    model.init_sims(replace=True)
    print("* successfully loaded the model !!!")
    return model

#model = load_model('/home/hclent/data/corpora/startrek/', 'startrek_model') ##print(model.syn0.shape) #(84196, 100)
#model = load_model('/home/hclent/data/18269575', 'coge_model') #(39605, 100) #100 features


def makeFeatureVec(words, model, num_features):
    # Function to average all of the word vectors in a given
    # paragraph
    #
    # Pre-initialize an empty numpy array (for speed)
    featureVec = np.zeros((num_features,),dtype="float32")
    #
    nwords = 0.
    #
    # Index2word is a list that contains the names of the words in
    # the model's vocabulary. Convert it to a set, for speed
    index2word_set = set(model.index2word)
    #
    # Loop over each word in the review and, if it is in the model's
    # vocabulary, add its feature vector to the total
    for word in words:
        if word in index2word_set:
            nwords = nwords + 1.
            featureVec = np.add(featureVec,model[word])
    #
    # Divide the result by the number of words to get the average
    featureVec = np.divide(featureVec,nwords)
    return featureVec



def getAvgFeatureVecs(num_docs, model, num_features):
    # Given a set of reviews (each one a list of words), calculate
    # the average feature vector for each one and return a 2D numpy array
    #
    # Initialize a counter
    counter = 0.
    #
    # Preallocate a 2D numpy array, for speed
    reviewFeatureVecs = np.zeros((len(num_docs),num_features),dtype="float32")
    #
    # Loop through the reviews
    for i in num_docs:
       #
       # Print a status message every 50th review
       if counter%50. == 0.:
           print ("Review %d of %d" % (counter, len(num_docs)))
       #
       # Call the function (defined above) that makes average feature vectors
       reviewFeatureVecs[counter] = makeFeatureVec(i, model, \
           num_features)
       #
       # Increment the counter
       counter = counter + 1.
    return reviewFeatureVecs



def run_word2vec(path_to_sentences, model_name, num_docs):
    sentences = MySentences(path_to_sentences)
    words = []
    for s in sentences:
        for word in s:
            words.append(word)
    model = load_model(path_to_sentences, model_name) #(39605, 100) #100 features
    featureVec = makeFeatureVec(words, model, 100)
    print("FEATURE VEC:")
    print(featureVec)
    num_docs = ['blah'] * num_docs
    reviewFeatureVecs = getAvgFeatureVecs(num_docs, model, 100)
    print("AVERAGE FEATURE VECTOR: ")
    print(reviewFeatureVecs)
    print("ONE ROW")
    average_vec = (reviewFeatureVecs[0]) #.tolist()
    print(type(average_vec))
    print(average_vec)
    print("SUM OF ONE ROW")
    print(np.sum(reviewFeatureVecs[0]))

    return average_vec

print("STAR TREK: ")
average_st = run_word2vec('/home/hclent/data/corpora/startrek/', 'startrek_model', 175)
print("-------------------------------------------------------------------------------")
print("COGE: ")
average_coge = run_word2vec('/home/hclent/data/18269575', 'coge_model', 165)

#print("DOT PRODUCT: ")
#normalized_st = LA.norm(average_st)
# normalized_coge = LA.norm(average_coge)

#print(np.dot(normalized_st, normalized_st))

