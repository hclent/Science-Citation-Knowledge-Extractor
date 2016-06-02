import nltk
from nltk import pos_tag, word_tokenize
import nltk.corpus
from nltk.corpus import stopwords
from collections import defaultdict
from collections import Counter

############ N GRAM FUNCTIONS #######################
def makeNgram(tagged_list, n, start_index):
	ngram = ""
	for nIdx in range(0, n):
		wordIdx = start_index + nIdx
		ngram = ngram + tagged_list[wordIdx] + " "
	ngram = ngram.strip()
	return ngram


def make_ALL_ngrams(tagged_list, n): #def ngram function #input the word_pos's and number
	ngram_list = [] #initialize a list
	for idx in range(len(tagged_list)-n):  #input word_pos and loop through words #stop before the end
		ngram = makeNgram(tagged_list, n, idx) #call ngram function
		ngram_list.append(ngram) #add our ngrams to the list
	return ngram_list


############## PRE PROCESSING ########################
def readMe(file):
	fcorpus = open(file, 'r')
	fcorpus = fcorpus.read()
	return fcorpus

#fcorpus = readMe("/Users/hclent/Desktop/webdev-biotool/flask/18088965_1.txt")
#print(fcorpus)


def formatCorpus(fcorpus): #as string
	#fcorpus = readMe(filehandle)
	fcorpus = fcorpus.replace("_", "underscore") #underscores will make problems for buildDict and are unimportant for my NER
	#make untagged corpus
	untagged_corpus = fcorpus.lower()
	untagged_corpus = untagged_corpus.split()
	stopwords = nltk.corpus.stopwords.words('english') #delete stopwords from untagged
	untagged_corpus = [w for w in untagged_corpus if w.lower() not in stopwords]

	fcorpus = fcorpus.lower()
	fcorpus = pos_tag(word_tokenize(fcorpus))

	#make tagged corpus and exclude stop words
	tag_corpus = []
	for wordAndTag in fcorpus:
		running_tag_corpus = []
		running_tag_corpus.append('_'.join(wordAndTag)) #word_tag

		#exclude stop words
		for wt in running_tag_corpus:
			(words, tags) = wt.split("_")
			if words not in stopwords:
				tag_corpus.append(wt)


	tag_corpus = ' '.join(tag_corpus)
	tag_corpus = tag_corpus.splitlines()


	return untagged_corpus, tag_corpus


# Word Boundaries
# def formBigrams(tag_corpus):
# 	bigramLines = []
# 	for oldLine in tag_corpus: #for the lines in the unigram set
# 	    newLines = "#_# " + oldLine # add #_# to the beginning for bigrams
# 	    bigramLines.append(newLines)
# 	return bigramLines


# # Word Boundaries
# def formTrigrams(tag_corpus):
# 	trigramLines = []
# 	for oldLine in tag_corpus:
# 	    newLines = "#_# #_# " + oldLine # #_# #_# word boundries for trigrams
# 	    trigramLines.append(newLines)
# 	return trigramLines


#Function for building ngram dictionaries, key:ngram, value:freq
def buildDict(n, format):
	wordsDict = defaultdict(lambda: 0)
	tagsDict = defaultdict(lambda: 0)
	taggedWordsDict = defaultdict(lambda: 0)
	sum = 0
	for sentence in format:
		taggedWords = sentence.split(" ")
		ngram_list = make_ALL_ngrams(taggedWords, n)
		for taggedNgram in ngram_list:
			if "_" in taggedNgram:
				word_gram = ""
				pos_gram = ""
				items = taggedNgram.split() #split the ngrams
				for item in items:
					(word, tag) = item.split("_") #split by word_tag #this will choke on URLs
					word_gram = word_gram + word + " "
					pos_gram = pos_gram + tag + " "
				word_gram = word_gram.strip()
				pos_gram = pos_gram.strip()
				wordsDict[word_gram] += 1
				tagsDict[pos_gram] += 1
				sum += 1
	return wordsDict, tagsDict


# unigrams, pos_unigrams = buildDict(1, corpus) #make unigram dictionary
# unigrams_dict = Counter(unigrams)
# pos_unigrams = Counter(pos_unigrams)

# bigrams, pos_bigrams = buildDict(2, bigramLines) #make bigram dictionary
# bigrams_dict = Counter(bigrams)
# pos_bigrams = Counter(pos_bigrams)

# trigrams, pos_trigrams = buildDict(3, trigramLines) #make trigram dictionary
# trigrams_dict = Counter(trigrams)
# pos_trigrams = Counter(pos_trigrams)


####################  MOST INFORMATIVE GRAMS ###############
# top_unigrams = unigrams_dict.most_common(100) #grab
# top_pos_unigrams = pos_unigrams.most_common(100)
# top_bigrams = bigrams_dict.most_common(100)
# top_pos_bigrams = pos_bigrams.most_common(100)
# top_trigrams = trigrams_dict.most_common(100)
# top_pos_trigrams = pos_trigrams.most_common(100)

# print ("The top 100 unigrams are: ")
# for unis in top_unigrams:
#     print (unis) #print the pos counts to it
# print ("---------------------------------------------")
# print ("The top 100 bigrams are: ")
# for bis in top_bigrams:
#     print (bis)
# print ("---------------------------------------------")
# print ("The top 100 trigrams are: ")
# for tris in top_trigrams:
#     print (tris)
