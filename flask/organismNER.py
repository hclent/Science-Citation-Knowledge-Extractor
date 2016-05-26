import nltk
from nltk import pos_tag, word_tokenize
from nltk.corpus import wordnet as wn
from nltk.corpus import stopwords
from collections import Counter
from collections import defaultdict
import re
from more_itertools import unique_everseen
import Ngrams as ngrams #myprogram
from collections import Counter
from nltk.stem.porter import PorterStemmer
from nltk.stem import WordNetLemmatizer


#Input: document(s)
#Output: NER dict for each doc


############ Feature Based Organism Named Entity Recognition ##############
#
#

named_entities = []

flags = ['genes', 'chromosomes', 'sequences', 'analysis', 'versus', 'searches', 'softwares'] #common false positives in latin NERs

#
#

#Finds two-word organism NE's with Latin based morphology #letter (dot) name
def matchOrganisms1(bigrams):
	organism_names = []
	for grams in bigrams:
		(w1, w2) = grams.split(" ")
		if w1[0].isalpha() and w1.endswith('.') and len(w1) == 2 and w2[-1].isalpha(): #also gets c.d. blah #BAD!
			organism = grams
			organism_names.append(organism) #['s. viridis', 'c. spinosa', 'c. gynandra']
			named_entities.append(organism)
	return organism_names

#
#

#Finds two-word organism NE's with Latin based morphology
def matchOrganisms2(bigrams):
	latin_names = []
	for grams in bigrams:
		(w1, w2) = grams.split(" ")
		p1 = re.compile('[a-z]{2,}[b-df-hj-np-tv-z]{1,}(i$|is$|ia$|ae$|um$|us$|es$|arum$)') #wont pick up ones that match in 'a'; too many false positives
		match1 = p1.search(w1)
		p2 = re.compile('[a-z]{2,}[b-df-hj-np-tv-z]{1,}(a$|i$|is$|ia$|ae$|um$|us$|es$|arum$)') #picks up the word 'genes'
		match2 = p2.search(w2)
		if match1 and match2:
			latin_names.append(grams)
			named_entities.append(grams)
	return (list(unique_everseen(latin_names)))

#
#

### Matches short latin names to long latin names for normalization
def normalizeNER(latinShort, latinLong, runningNErlist):  #Right now this is ONE dictionary per document
	runningDict = buildDict(runningNErlist)
	dictLatin = buildDict(latinLong)

	for short_names in latinShort: #list
		(w1, w2) = short_names.split(" ") #split by bigram
		for long_names in latinLong:
			(x1, x2) = long_names.split(" ")
			if w2 == x2:
				#if the second word in each latin name are the same
				#transfer counts so that we have one value per organism
				#print("MATCHES:")
				#print(str(short_names) + " === " + str(long_names)) ##s. italica === setaria italica
				dictLatin[long_names] += runningDict[short_names] #update dict #add short name count to long name count
				#Nothing's been changed in runningDict....
				#
			#else:
				# if there is only the short name for something, we need to keep it #add it to dict
				#dictLatin[short_names] += runningDict[short_names] #PROBLEM: THIS IS WONKY
	#print("NORMALIZED DICT")
	#print(dictLatin)
	#print("LENGTH OF DICT: ")
	#print(len(dictLatin))


#
#

### Finds one-word organism NE's based on definitions and similarity score with WordNet gazetteer
def wordNetNER(document):
	plant_sns = (wn.synsets('plant', pos="n"))
	plant = plant_sns[1] #(botany) a living organism lacking the power of locomotion #hardcoded

	wordnet_names = []
	#wordnet_lemmatizer = WordNetLemmatizer   ##Lematizer doesn't work....
	for word in document:
		#word = wordnet_lemmatizer.lemmatize(word) ##Lematizer doesn't work...
		mySynsets = wn.synsets(word, pos="n")

		i = 0
		for i in range(0, 3):
			try:
				given_word = mySynsets[i] #tries first 3 synsets
				definition = (given_word.definition())
				p1 = re.compile('plant(s?)\s')
				p2 = re.compile('organism(s?)\s')
				p3 = re.compile('animal(s?)\s')
				match1 = p1.search(definition)
				match2 = p2.search(definition)
				match3 = p3.search(definition)

				if match1 or match2 or match3:  #if the given word has "plants" or "animals" in the def, check to see how similar it is to "plant"
					similarity_score = (given_word.path_similarity(plant)) #check similarity score
					if similarity_score >= 0.2:
						#print(similarity_score)
						#print ("The words: "+(str(given_word)) + "  has a sim score of:  " +str(similarity_score))
						wordnet_names.append(word)
						named_entities.append(word)
			#hypernym = given_word.hypernyms() #hypernym is list #synset 'organism' exists #can't search in the hypernyms....hmm...
				i += 1
			except IndexError:
				pass
	wordnet_ner = (list(unique_everseen(wordnet_names)))
	return wordnet_ner

#
#


#### For counts
def buildDict(list):
	wordsDict = defaultdict(lambda: 0)
	sum = 0
	for ner in list:
		wordsDict[ner] +=1
		sum += 1
	return wordsDict

#
#

######## LOADING AND PRE PROCESSING TEXT ##########
#prefix = "/PycharmProjects/untitled/Crawling/"

def loadDocuments(filenamePrefix, maxNum):
	#print(" * load messages: Started .... ")
	for i in range(1, maxNum+1):
		#print("MESSAGE: " +str(i))

		filename = filenamePrefix +'_'+ str(i) + ".txt"
		fcorpus = ngrams.readMe(filename) #Function from Ngrams.py
		untagged_corpus, tagged_corpus = ngrams.formatCorpus(fcorpus) #No stop words


		bigramLines = [] #No stop words
		for oldLine in tagged_corpus:
			newLines = "#_# " + oldLine
			bigramLines.append(newLines)
		bigrams, pos_bigrams = ngrams.buildDict(2, bigramLines)
		bigrams_dict = Counter(bigrams)

		#print(" * load messages: complete !!!")
		#print(" * Organism Named Entity Recognition: Started .... ")

		latin_short = matchOrganisms1(bigrams_dict)
		print("LATIN SHORT:")
		print(latin_short)

		latin_long = matchOrganisms2(bigrams_dict)
		print("LATIN LONG:")
		print(latin_long)

		normalizeNER(latin_short, latin_long, named_entities)

		wordnet_names = wordNetNER(untagged_corpus)
		print("WORDNET: ")
		print(wordnet_names)

		#print(" * Organism NER: complete!!!" + "\n")

	#nerDict = Counter(buildDict(named_entities))
	#print(len(nerDict))
	#return nerDict
	return latin_short, latin_long, wordnet_names

#loadDocuments(prefix, 10)

#
#

######## NOTES #####
# senses = (wn.synsets('plant', pos="n"))
# s1 = senses[0]
# s2 = senses[1]
# s3 = senses[2]
# s4 = senses[3]
# d1 = (s1.definition()) #buildings for carrying on industrial labor
# d2 = (s2.definition()) #(botany) a living organism lacking the power of locomotion
# d3 = (s3.definition()) #an actor situated in the audience whose acting is rehearsed but seems spontaneous to the audience
# d4 = (s4.definition()) #something planted secretly for discovery by another


#TO DO:
# 1. fix nerNormalize problems
# 2. Handling typos within the academic journals
# 3. Do lemmatization
# 4. See how hypernymns perform as a feature
# 5. modify output for a list
# Should have:
#[[Author, titles, urls, journal, [NERs]],[Doc2 info],[Doc3 info]...]
#
# 6. Rank NER's based on freq to use as feature in classifier
