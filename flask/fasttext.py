#### Playspace for clustering on fasttext embeddings
from database_management import retrieveAllPmids
from multi_preprocess import retrieveBioDocs
from processors import *
import sys, os.path, time
import numpy as np
from collections import defaultdict, Counter
import gensim, os, codecs
from gensim.models import Word2Vec
from sklearn.cluster import KMeans
import math, random



'''
Step 1:
    Get lemmas and tags from biodocs
    [[doc1 words], [doc2 words]] & [[doc1 tags],[doc2 tags]]
Step 2:
    NP tokenization
    [[doc1 tokens], [doc2 tokens]]
Step 3:
    flatmap tokens for training w2v
    [all tokens ever]
    train w2v on all tokens
Step 4:
    Try these:
    Get the top X% NPs,
    Rank NPs with TF-IDF,
Step 5:
    concatenate word embeddings to make doc embeddings for each doc
Step 6:
    cluster
'''


def flatten(listOfLists):
    return list(chain.from_iterable(listOfLists))

def get_words_tags(pmid_list):
    words = []
    tags = []

    for pmid in pmid_list:
        biodocs = retrieveBioDocs(pmid)
        print(len(biodocs))

        for doc in biodocs:
            doc_words = []
            doc_tags = []
            with open(doc) as jf:
                data = Document.load_from_JSON(json.load(jf))
                num_sentences = data.size
                for i in range(0, num_sentences):
                    s = data.sentences[i]
                    s_words = s.lemmas
                    s_tags = s.tags
                    doc_words.append(s_words)
                    doc_tags.append(s_tags)
            doc_words = flatten(doc_words)
            doc_tags = flatten(doc_tags)
            words.append(doc_words)
            tags.append(doc_tags)

    #words = flatten(words)
    #tags = flatten(tags)
    return words, tags


def transform_text(words, tags):
    transformed_tokens = []
    np_stack = []
    join_char = "_"
    for i, w in enumerate(words):
        if tags[i].startswith("NN"):
            np_stack.append(w)
        else:
            if len(np_stack) > 0:
                np = join_char.join(w.lower() for w in np_stack)
                transformed_tokens.append(np)
                # reset stack
                np_stack = []
            # append current token
            transformed_tokens.append(w)
    return transformed_tokens


def chooseTopNPs(transformed_tokens):
    npDict = Counter(defaultdict(lambda: 0))
    for t in transformed_tokens:
        if "_" in t:
            npDict[t] += 1
    return npDict


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

    # if you need to add sentences:
    #more_sentences = MySentences('/home/hclent/data/18952863')
    #model.train(more_sentences)

    model.save(str(path_to_sentences+'/'+model_name))
    print("* successfully saved the model !!!")


def load_model(file):
    print("* loading the model ... ")
    if ".vec" not in file:
        model = Word2Vec.load(str(file))
    if ".vec" in file:
        model = Word2Vec.load_word2vec_format(file, binary=False)  # C binary format
    model.init_sims(replace=True)
    print("* successfully loaded the model !!!")
    return model



def makeFeatures(words, tags, top_nps, model):
    feature_vecs = []
    size = model.syn0.shape
    num_features = size[1]
    print(top_nps)
    i = 0
    for doc in words:
        doc_vec = []
        transformed = transform_text(doc, tags[i])
        #print(transformed)
        for noun in top_nps: #top_ns is a counter
            #print(noun)

            if noun[0] in transformed:
                np_vec = model[noun[0]]
                #print(np_vec)
                doc_vec.append(np_vec)
            else:
                v_list = [0.0] * int(num_features)
                empty_vec = np.array(v_list)
                #print(empty_vec)
                doc_vec.append(empty_vec)
                #append a vec of length num_features
        doc_vec = np.array(flatten(doc_vec))
        feature_vecs.append(doc_vec)

        i += 1

    print(len(feature_vecs))
    print(type(feature_vecs[0]))
    print(feature_vecs[0])
    print(len(feature_vecs[0]))
    return feature_vecs


# make list of embeddings to cluster [[embedding 1], [embedding 2], ... ]
def getNPvecs(top_nps, model):
    all_vecs = []
    for nps in top_nps:
        np_vec = model[nps[0]]
        vec = np_vec.tolist()
        all_vecs.append(vec)
    matrix = np.array(all_vecs)
    return matrix






# pmid_list = retrieveAllPmids()
# words, tags = get_words_tags(pmid_list) #list of words/tags per doc
# t1 = time.time()
# flat_words = flatten(words)
# flat_tags = flatten(tags)
# print("with flatmap " + str(time.time() - t1))
# transformed_sentence = transform_text(flat_words, flat_tags)
# npDict = chooseTopNPs(transformed_sentence)
# top_nps = list(npDict.most_common(300))
# print(top_nps)
top_nps = [('gene_expression', 3448), ('t_\\_t', 3421), ('c._jejunus', 2433), ('amino_acid', 2247), ('climate_change', 2101), ('expression_level', 2021), ('cell_cycle', 2006), ('stem_cell', 1973), ('scale_bar', 1943), ('transcription_factor', 1912), ('cell_line', 1789), ('cell_death', 1784), ('\\_t', 1703), ('b._napus', 1680), ('gene_family', 1680), ('t_cell', 1612), ('cell_type', 1581), ('time_point', 1545), ('a._thaliana', 1494), ('plant_species', 1454), ('breast_cancer', 1417), ('dna_methylation', 1374), ('cancer_cell', 1335), ('land_plant', 1326), ('error_bar', 1321), ('expression_pattern', 1305), ('dna_damage', 1290), ('signaling_pathway', 1271), ('room_temperature', 1226), ('cell_proliferation', 1191), ('table_s1', 1168), ('candidate_gene', 1121), ('target_gene', 1097), ('_', 1067), ('animal_model', 1042), ('western_blot', 1016), ('motor_neuron', 1007), ('mouse_model', 988), ('tumor_cell', 987), ('genome_sequence', 968), ('datum_set', 940), ('cell_division', 935), ('plant_genome', 928), ('protein_level', 923), ('gene_duplication', 888), ('growth_factor', 885), ('arabidopsis_thaliana', 874), ('b_cell', 865), ('s_phase', 857), ('bone_marrow', 844), ('\\_usepackage', 840), ('stress_granule', 836), ('cell_wall', 830), ('dna_replication', 801), ('plasma_membrane', 801), ('b._rapa', 787), ('zebra_finch', 766), ('western_blotting', 747), ('table_s2', 736), ('gene_tree', 735), ('western_blot_analysis', 710), ('protein_sequence', 706), ('progenitor_cell', 703), ('material_online', 702), ('risk_factor', 691), ('cell_growth', 681), ('amino_acid_sequence', 659), ('%_ci', 658), ('stress_response', 645), ('pcr_product', 642), ('control_group', 639), ('flowering_plant', 639), ('expression_profile', 631), ('flow_cytometry', 631), ('cell_cycle_progression', 609), ('gene_pair', 604), ('genome_duplication', 599), ('hh_signaling', 594), ('p._patens', 587), ('sequence_alignment', 586), ('brain_region', 577), ('reference_genome', 574), ('sequence_similarity', 570), ('duplication_event', 558), ('figure_supplement', 554), ('\\_u2009', 554), ('protein_kinase', 534), ('genome_size', 532), ('p_value', 531), ('mm_nacl', 530), ('protein_expression', 525), ('table_s3', 521), ('species_tree', 514), ('mrna_level', 512), ('rna_binding_protein', 500), ('crystal_structure', 500), ('linkage_group', 499), ('binding_site', 499), ('copy_number', 488), ('seed_plant', 485)]

#
#
#
#
#
# print("making transformed s")
# transformed_s = (' '.join(transformed_sentence))
# print(transformed_s[:1000])
# suffix = "/home/hclent/data/nns/5k"
# completeName = os.path.join(suffix, ('5210_nns.txt'))  #pmcid.txt #save to suffix path
# sys.stdout = open(completeName, "w")
# print(transformed_s)

# create_model('/home/hclent/data/nns/5k', '5210_nns')
def get_hashing(data):
  t0 = time.time()
  hasher = HashingVectorizer() #l2 projected on the euclidean unit sphere
  hX = hasher.fit_transform(data)
  return hX, hasher


#Input: High dimensional (sparse) matrix
#Output: Clusters
# labels = km.labels_
# centroids = km.cluster_centers_
def do_kemeans(sparse_matrix, k_clusters):
    t0 = time.time()
    num_clusters = int(k_clusters)
    km = KMeans(init='k-means++', n_clusters=num_clusters)
    km.fit(sparse_matrix)
    clusters = km.labels_.tolist()
    return clusters


model = load_model('/home/hclent/tmp/fastText/ftmodel.vec')
#model = load_model("/home/hclent/data/nns/5k/5210_nns")
print(model.syn0.shape)
matrix = getNPvecs(top_nps, model)
kmeans = KMeans(n_clusters=10, random_state=2).fit(matrix)


res = list(zip(kmeans.labels_, top_nps))

t1 = []
t2 = []
t3 = []
t4 = []
t5 = []
t6 = []
t7 = []
t8 = []
t9 = []
t10 = []
for t in res:
    if t[0] == 0:
        t1.append(t[1])
    if t[0] == 1:
        t2.append(t[1])
    if t[0] == 2:
        t3.append(t[1])
    if t[0] == 3:
        t4.append(t[1])
    if t[0] == 4:
        t5.append(t[1])
    if t[0] == 5:
        t6.append(t[1])
    if t[0] == 6:
        t7.append(t[1])
    if t[0] == 7:
        t8.append(t[1])
    if t[0] == 8:
        t9.append(t[1])
    if t[0] == 9:
        t10.append(t[1])
print("t1: ")
print(t1)
print("t2: ")
print(t2)
print("t3: ")
print(t3)
print("t4: ")
print(t4)
print("t5: ")
print(t5)
print("t6: ")
print(t6)
print("t7: ")
print(t7)
print("t8: ")
print(t8)
print("t9: ")
print(t9)
print("t10: ")
print(t10)




# fv = makeFeatures(words, tags, top_nps, model)
# print(fv)
#
#
# hX, hasher = get_hashing(fv) #returns hX and hasher
# clusters = do_kemeans(hX, 10) #returns clusters
#
# rank = zip(list(top_nps, clusters))
# println(rank)
