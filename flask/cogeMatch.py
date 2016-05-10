import nltk
from nltk import pos_tag, word_tokenize
from collections import Counter
from collections import defaultdict
import re
import Ngrams as ngrams #myprogram

## What CoGe tools are being used the most?
## Input MainCrawler's text files
## Output's a dictionary with counts
#### ## add in for 'Document doesn't say' #####

toolsDict = defaultdict(lambda: 0)

def findTools(unigrams, bigrams):

	for grams in unigrams:
		synmap = re.compile('synmap')
		synfind = re.compile('synfind')
		gevo = re.compile('gevo')
		featview = re.compile('featview')

		match1 = synmap.search(grams)
		match2 = synfind.search(grams)
		match3 = gevo.search(grams)
		match4 = featview.search(grams)

		if match1:
			toolsDict['synmap'] += 1
		else:
			toolsDict['synmap'] == 0
		if match2:
			toolsDict['synfind'] += 1
		else:
			toolsDict['synfind'] == 0
		if match3:
			toolsDict['gevo'] += 1
		else:
			toolsDict['gevo'] == 0
		if match4:
			toolsDict['featview'] += 1
		else:
			toolsDict['featview'] == 0

	for two_grams in bigrams:
		(w1, w2) = two_grams.split(" ")
		organism = re.compile('organism')
		view = re.compile('view')

		match5 = organism.search(w1)
		match6 = view.search(w2)
		if match5 and match6:
			toolsDict['organism view'] += 1
		else:
			toolsDict['organism view'] == 0

		coge = re.compile('coge')
		blast = re.compile('blast')

		match7 = coge.search(w1)
		match8 = blast.search(w2)
		if match7 and match8:
			toolsDict['coge blast'] +=1
		else:
			toolsDict['coge blast'] == 0


	return toolsDict



prefix = "/PycharmProjects/untitled/Crawling/"

def loadDocuments(filenamePrefix, maxNum):
	print(" * load messages: Started .... ")
	for i in range(1, maxNum+1):
		print("MESSAGE: " +str(i))

		filename = filenamePrefix + str(i) + ".txt"
		fcorpus = ngrams.readMe(filename) #Function from Ngrams.py

		untagged_corpus, corpus = ngrams.formatCorpus(fcorpus)
		tagged_corpus = corpus.splitlines()

		unigrams, pos_unigrams = ngrams.buildDict(1, tagged_corpus) #Function from Ngrams.py

		bigramLines = []
		for oldLine in tagged_corpus:
			newLines = "#_# " + oldLine # add #_# to the beginning for bigrams # '#' for word boundaries
			bigramLines.append(newLines)

		bigrams, pos_bigrams = ngrams.buildDict(2, bigramLines) #make bigram dictionary

		findTools(unigrams, bigrams)

loadDocuments(prefix, 145)

print(toolsDict)
