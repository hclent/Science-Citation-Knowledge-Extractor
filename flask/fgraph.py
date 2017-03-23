from gensim.models import Word2Vec
import time
import numpy as np

t1 = [('gene_expression', 3448), ('expression_level', 2021), ('cell_line', 1789), ('time_point', 1545), ('expression_pattern', 1305)]
t2 = [('transcription_factor', 1912), ('signaling_pathway', 1271), ('target_gene', 1097), ('growth_factor', 885), ('stress_response', 645)]
t3 = [('climate_change', 2101), ('gene_family', 1680), ('plant_genome', 928), ('gene_duplication', 888), ('gene_tree', 735)]

topics = [t1, t2, t3]

#load fasttext file
def load_model(file):
    print("* loading the model ... ")
    t0 = time.time()
    if ".vec" not in file:
        model = Word2Vec.load(str(file))
    if ".vec" in file:
        model = Word2Vec.load_word2vec_format(file, binary=False)  # C binary format
    model.init_sims(replace=True)
    print("* successfully loaded the model: done in %0.3fs." % (time.time() - t0))
    print(model.syn0.shape)
    return model


#make similarity matrix to decide what is/isn't connected in force dir graph
#words in same topic are connected (1)
#words in different topics are connected if cos(w1,w2) > some x
def make_matrix(topics, model):
    matrix = []
    for t in topics:

        for word in t:
            temp_row = []
            for i in topics:
                if t == i :
                    for word2 in i:
                        if word == word2:
                            temp_row.append(word[1])
                        if word != word2:
                            temp_row.append(1.0)
                if t != i:
                    for word2 in i:
                        print("cos(" + str(word[0]) + ", " + str(word2[0]) + ")")
                        np_float = model.similarity(word[0], word2[0])
                        w_score = np_float.item()
                        temp_row.append(w_score)
            matrix.append(temp_row)
    return matrix


model = load_model('/home/hclent/tmp/fastText/ftmodel.vec')
matrix = make_matrix(topics, model)

words_list = []

print("TOPICS")
for t in topics:
    for words in t:
        print(words[0])
        words_list.append(words[0])

print((', ').join(words_list))