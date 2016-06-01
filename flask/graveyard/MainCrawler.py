import requests
from bs4 import BeautifulSoup
import re
import sys
import pickle
import socket
import socks

## Run in Python3 for handling unicode!
## Input: PMC 'cited by' page
## Output: Retrieves all article titles, authors, url, journal of publication, and document text

## May 31, 2016: This code has become outdated and is no longer used to do information retrieval
## May 31, 2016: See Entrez_IR.py for new information retrieval system

def pmc_spider(max_pages, pmid): #Main spider
	start = 1

	titles_list = []
	url_list = []
	url_keys = []
	authors_list = []

	while start <= max_pages:
		url = 'http://www.ncbi.nlm.nih.gov/pmc/articles/pmid/'+str(pmid)+'/citedby/?page='+str(start)

		req = requests.get(url)
		plain_text = req.text
		soup = BeautifulSoup(plain_text, "lxml")

		for items in soup.findAll('div', {'class': 'title'}):
			title = items.get_text()
			titles_list.append(title)

			for link in items.findAll('a'):
				urlkey = link.get('href')
				url_keys.append(urlkey)   #url = base + key
				url =  "http://www.ncbi.nlm.nih.gov"+str(urlkey)
				url_list.append(url)

		for citation in soup.findAll('div', {'class': 'desc'}):
			author = citation.string
			authors_list.append(author)

		start += 1
	return titles_list, url_list, authors_list


titles_list, url_list, authors_list = pmc_spider(1, '18269575') #look at first n pages



def get_text(list_paper_urls, pmid): #Get publication text
	print ("* Retrieving papers .... this will take a while .... ")

	i = 1
	academic_journals = []

	for paper_url in list_paper_urls:
		#print("* accessing paper from......")
		#print(paper_url)
		main_text = []

		#print("* getting requests .... ")
		req = requests.get(paper_url)

		#print("* request obtained .... ")
		plain_text = req.text
		soup = BeautifulSoup(plain_text, 'lxml')

		# Get journal of publication #Slow
		for journals in soup.findAll('a', class_ = 'navlink', href = re.compile("journal")):
			journal_name = journals.get_text()
			not_a_journal = re.compile('Journal List') #exclude this unwanted text
			discard = not_a_journal.search(journal_name)
			if not discard:
				academic_journals.append(journal_name)

		#Get main text
		for words in soup.findAll('p', id = re.compile("__p"), class_ = re.compile('first')):
			document = words.get_text()
			main_text.append(document)
		main_text = ' '.join(main_text)
		print(main_text)

		#print to .txt file
		#sys.stdout = open(str(pmid)+'_'+str(i)+'.txt', "w")   ## change "w" to "r+" if file exists
		#print (main_text)

		i += 1
	return academic_journals



#print(len(url_count))
#print(len(main_info)) 
#if these numbers don't match, the bot has been blocked on some pages, rerun with proxies


############ GRAVEYARD #############
#keywords = []
		# for abstracts in soup.findAll('span', {'class': 'kwd-text'}): #keywords #will use as feature
		# 	keyws = abstracts.get_text()
		# 	keywords.append(keyws)
		# 	print(keywords)
# I can't recall where I'm getting these from on a page, and they don't seem to match for the papers...
