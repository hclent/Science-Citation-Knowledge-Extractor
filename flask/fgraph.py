from gensim.models import Word2Vec
import time
import re

t1 = [('gene_expression', 3448), ('expression_level', 2021), ('cell_line', 1789), ('time_point', 1545), ('expression_pattern', 1305), ('western_blot', 1016), ('protein_level', 923), ('western_blotting', 747), ('western_blot_analysis', 710), ('pcr_product', 642), ('expression_profile', 631), ('protein_expression', 525), ('mrna_level', 512)]
t2 = [('transcription_factor', 1912), ('signaling_pathway', 1271), ('target_gene', 1097), ('growth_factor', 885), ('stress_response', 645), ('protein_kinase', 534), ('rna_binding_protein', 500)]
t3 = [('climate_change', 2101), ('gene_family', 1680), ('genome_sequence', 968), ('datum_set', 940), ('plant_genome', 928), ('gene_duplication', 888), ('gene_tree', 735), ('protein_sequence', 706), ('material_online', 702), ('gene_pair', 604), ('genome_duplication', 599), ('reference_genome', 574), ('sequence_similarity', 570), ('duplication_event', 558), ('genome_size', 532), ('species_tree', 514), ('copy_number', 488)]
t4 = [('c._jejunus', 2433), ('cell_cycle', 2006), ('cell_death', 1784), ('dna_damage', 1290), ('room_temperature', 1226), ('cell_proliferation', 1191), ('cell_division', 935), ('s_phase', 857), ('stress_granule', 836), ('dna_replication', 801), ('cell_growth', 681), ('flow_cytometry', 631), ('cell_cycle_progression', 609)]
t5 = [('t_\\_t', 3421), ('scale_bar', 1943), ('\\_t', 1703), ('error_bar', 1321), ('_', 1067), ('\\_usepackage', 840), ('%_ci', 658), ('\\_u2009', 554), ('p_value', 531), ('mm_nacl', 530)]
t6 = [('stem_cell', 1973), ('t_cell', 1612), ('cell_type', 1581), ('breast_cancer', 1417), ('dna_methylation', 1374), ('cancer_cell', 1335), ('tumor_cell', 987), ('b_cell', 865), ('bone_marrow', 844), ('progenitor_cell', 703), ('control_group', 639), ('hh_signaling', 594)]
t7 = [('amino_acid', 2247), ('plasma_membrane', 801), ('amino_acid_sequence', 659), ('sequence_alignment', 586), ('figure_supplement', 554), ('crystal_structure', 500), ('binding_site', 499)]
t8 = [('b._napus', 1680), ('table_s1', 1168), ('candidate_gene', 1121), ('b._rapa', 787), ('table_s2', 736), ('table_s3', 521), ('linkage_group', 499)]
t9 = [('animal_model', 1042), ('motor_neuron', 1007), ('mouse_model', 988), ('zebra_finch', 766), ('risk_factor', 691), ('brain_region', 577)]
t10 =[('a._thaliana', 1494), ('plant_species', 1454), ('land_plant', 1326), ('arabidopsis_thaliana', 874), ('cell_wall', 830), ('flowering_plant', 639), ('p._patens', 587), ('seed_plant', 485)]


def make_lut(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10):
    groupDict = {}
    for word_count in t1:
        groupDict[word_count[0]] = '1'
    for word_count in t2:
        groupDict[word_count[0]] = '2'
    for word_count in t3:
        groupDict[word_count[0]] = '3'
    for word_count in t4:
        groupDict[word_count[0]] = '4'
    for word_count in t5:
        groupDict[word_count[0]] = '5'
    for word_count in t6:
        groupDict[word_count[0]] = '6'
    for word_count in t7:
        groupDict[word_count[0]] = '7'
    for word_count in t8:
        groupDict[word_count[0]] = '8'
    for word_count in t9:
        groupDict[word_count[0]] = '9'
    for word_count in t10:
        groupDict[word_count[0]] = '10'
    return groupDict


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
#words are assigned a group id for the topic they fall in
def make_matrix(topics, model, groupDict, words_list):
    val_matrix = []

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
                        #print("cos(" + str(word[0]) + ", " + str(word2[0]) + ")")
                        np_float = model.similarity(word[0], word2[0])
                        w_score = np_float.item()
                        temp_row.append(w_score)
            val_matrix.append(temp_row)

            #group val is last val on each row
            group_num = int(groupDict.get(word[0]))
            temp_row.append(group_num)

    #add final rows column
    last_row_vals = [ int(groupDict.get(word)) for word in words_list ]
    val_matrix.append(last_row_vals)

    return val_matrix


model = load_model('/home/hclent/tmp/fastText/ftmodel.vec')
topics = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10]
groupDict = make_lut(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10)
words_list = []
for t in topics:
    for words in t:
        #print(words[0])
        words_list.append(words[0])


val_matrix = make_matrix(topics, model, groupDict, words_list)


def make_csv(val_matrix, words_list):
    top_row = str(str(",")+(",".join(words_list)) + str(",group"))
    print(top_row)
    i = 0
    rows = words_list + ["group"]
    for label in rows: #label is str
        values = val_matrix[i] #list
        line = str( label + "," + str(values))
        line = re.sub( "\s+" , "", line)
        line = re.sub("\]|\[", "", line)
        print(line)
        i+=1


make_csv(val_matrix, words_list)