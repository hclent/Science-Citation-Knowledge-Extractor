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
main_info = list(zip(titles_list, authors_list))
#url_count = (len(url_list))

#Finding and Comparing Syntenic Regions among Arabidopsis
#pmid = 18952863

#How to usefully compare homologous plant genes and chromosomes as DNA sequences.
#pmid = 18269575


####### PICKLE OBJECTS ############

# pickle URL list
#with open('maincrawler_urls_list.pickle', 'wb') as f:
#	pickle.dump(url_list, f) #


# pickle authors list
#with open('maincrawler_authors.pickle', 'wb') as f:
#	pickle.dump(authors_list, f)


# pickle titles 
#with open('maincrawler_titles.pickle', 'wb') as f:
#	pickle.dump(titles_list, f)


pickle_in = open('maincrawler_urls_list.pickle', 'rb')
url_pickle = pickle.load(pickle_in)


pickle_in2 = open('maincrawler_authors.pickle', 'rb')
authors_pickle = pickle.load(pickle_in2)


pickle_in3 = open('maincrawler_titles.pickle', 'rb')
titles_pickle = pickle.load(pickle_in3)
#########################################


#Now takes pmid also for naming

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


#journals_list = get_text(url_pickle)


####### PICKLE OBJECTS ############
#with open('maincrawler_journals.pickle', 'wb') as f:
#	pickle.dump(journals_list, f)

pickle_in4 = open('maincrawler_journals.pickle', 'rb')
journals_pickle = pickle.load(pickle_in4)
####################################

main_pickle_info = list(zip(titles_pickle, authors_pickle, journals_pickle, url_pickle))

#print(len(url_count))
#print(len(main_info)) 
#if these numbers don't match, the bot has been blocked on some pages. 
#rerun with proxy


############ GRAVEYARD #############
#keywords = []
		# for abstracts in soup.findAll('span', {'class': 'kwd-text'}): #keywords #will use as feature
		# 	keyws = abstracts.get_text()
		# 	keywords.append(keyws)
		# 	print(keywords)
# I can't recall where I'm getting these from on a page, and they don't seem to match for the papers...
