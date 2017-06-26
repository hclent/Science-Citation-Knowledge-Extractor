from scipy import spatial
import numpy as np


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

