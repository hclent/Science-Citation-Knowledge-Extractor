import gensim, logging, os, codecs
from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence
from nltk.corpus import brown


logging.basicConfig(filename='.2vec.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
logging.info('Started')


#read startrek texts
#tokenize them into sents somehow
#have sentences = [['word', 'words'],[]]
#feed those into model

class MySentences(object):
    def __init__(self, dirname):
        self.dirname = dirname

    def __iter__(self):
        for fname in os.listdir(self.dirname):
            if ".txt" in fname:
                for line in codecs.open(os.path.join(self.dirname, fname), "r", encoding='utf-8', errors='ignore'):
                    if len(line.split()) != 0: #make sure its not empty
                        yield line.split()


def create_model(path_to_sentences, model_name):
    print("* retrieving sentences ... ")
    sentences = MySentences(path_to_sentences) # a memory-friendly iterator
    print("* successfully retrieved sentences !!!")
    print("* making the model ... ")
    model = gensim.models.Word2Vec(sentences, min_count=1)
    print("* successfully made the model !!!")
    print("* saving the model ... ")
    model.save(str(path_to_sentences+model_name))
    print("* successfully saved the model !!!")

#create_model('/home/hclent/data/corpora/startrek/', 'startrek_model')

def load_model(path_to_sentences, model_name):
    print("* loading the model ... ")
    model = Word2Vec.load(str(path_to_sentences+model_name))
    model.init_sims(replace=True)
    print("* successfully loaded the model !!!")
    return model

model = load_model('/home/hclent/data/corpora/startrek/', 'startrek_model')


print(model.most_similar(['Borg'], topn=3))
print(model.doesnt_match("Picard Worf Riker Borg".split()))


more_examples = ["Picard Borg Data",  "Enterprise Picard cube"]
for example in more_examples:
    a, b, x = example.split()
    predicted = model.most_similar([x, b], [a])[0][0]
    print("'%s' is to '%s' as '%s' is to '%s'" % (a, b, x, predicted))

# print("Woman is to Girl, as Man is to .... ")
# print(model.most_similar(positive=['woman', 'girl'], negative=['man'], topn=3)) #precomputing L2-norms of word weight vectors




