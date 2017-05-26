import time, datetime
from collections import defaultdict
from sqlalchemy import create_engine, MetaData, Table, select
#from sqlalchemy.sql import select
# from sqlalchemy.orm import scoped_session, sessionmaker
# from sqlalchemy.ext.declarative import declartive_base


#May 26, 2017: accessing db through sqlite3 library depreciated



'''
Tables in pmids_info.db
1) inputPapers: stores info about the input pmids and their citations
2) citations: stores information about pmcids that cite input pmids
3) queries: stores information about a particular input query e.g. paper1+paper2
4) annotations: stores information about the annotated pmcids' lemmas and named entities categories

'''


#Connect to db through SqlAlchemy:
engine = create_engine('sqlite:///pmids_info.db')

def connection():
	conn = engine.connect()
	return conn

def load_tables():
	metadata = MetaData(bind=engine) #init metadata. will be empty
	metadata.reflect(engine) #retrieve db info for metadata (tables, columns, types)
	inputPapers = Table('inputPapers', metadata)
	citations = Table('citations', metadata)
	queries = Table('queries', metadata)
	annotations = Table('annotations', metadata)
	return inputPapers, citations, queries, annotations

inputPapers, citations, queries, annotations = load_tables()



#################### SUPPORT FUNCTIONS FOR inputPapers TABLE ######################
#Input: pmid
#Output: apa citation of *THAT* pmid
#Updated to use sqlalchemy
def db_inputPapers_retrieval(user_input):
	result = engine.execute('select title, author, journal, pubdate, url from inputPapers where pmid = :0', [user_input])
	for row in result:
		title = row['title']
		author = row['author']
		journal = row['journal']
		pubdate = row['pubdate']
		url = row['url']
		apa = str(author+' ('+pubdate+'). '+title+'. '+journal+'. Retrieved from '+url)
		return apa


#Input: pmid
#Output: number of citations
#Updated to use sqlalchemy
def db_input_citations_count(user_input):
	result = engine.execute('select num_citations from inputPapers where pmid = :0', [user_input])
	for row in result:
		num_citations = row['num_citations']
		return num_citations


# Input: pmid that is cited
# Output: dicts needed for statistics tab
# Updated to use sqlalchemy
def db_statistics(user_input):
	pmidDict = defaultdict(int)
	pmcDict = defaultdict(list)

	num_citations = db_input_citations_count(user_input)
	pmidDict[user_input] = num_citations

	result = engine.execute('select pmcid, abstract, whole_article, sents, tokens from citations where citesPmid = :0', [user_input])
	for row in result:
		pmcid = row['pmcid']
		abstract = row['abstract']
		whole = row['whole_article']
		sents = row['sents']
		tokens = row['tokens']
		pmcDict[pmcid] = [abstract, whole, sents, tokens]

	return pmidDict, pmcDict



# After you've scraped the input_paper, write the information about abstract, article, and self_pmcid to db
# from function getSelfText(user_input)
# puts pmcid, abstract check, and article check into db from
# Updated to sqlalchemy syntax
def updateInputPapers(user_input, self_pmcid, abstract, article):
	up = inputPapers.update().\
		where(inputPapers.c.pmid == user_input).\
		values(dict(abstract=abstract, whole_article=article, pmcid=self_pmcid))
	conn = connection()
	conn.execute(up)



# convert pmid2pmcid
# updated to sqlAlchemy
def pmid2pmcid(user_input):
	result = engine.execute('select pmcid from inputPapers where pmid = :0',[user_input])
	for row in result:
		pmcid = row["pmcid"]  # return first thing in tuple ('2836516',)
		return pmcid
	# will return NoneType if its not there apparently :)



######################### SUPPORT FUNCTIONS FOR citations TABLE ###########################
#Input: pmid
#Output: apa citations for citing pmCids as hyperlinks
#Updated to sqlAlchemy
def db_citations_hyperlink_retrieval(user_input):
	result = engine.execute('select title, author, journal, pubdate, url from citations where citesPmid = :0', [user_input])
	apa_citations = []
	for row in result:
		title = row["title"]
		author = row["author"]
		journal = row["journal"]
		pubdate = row["pubdate"]
		url = row["url"]
		apa = str(author+' ('+pubdate+'). '+title+'. '+journal+'. Retrieved from '+url)
		href_label = str('<a href="'+url+'">'+str(apa)+'</a>')
		apa_citations.append(href_label)
	return apa_citations


#Input: pmid
#Output: list of apa citations of pmc-ids citing that pmid
#Updated to sqlAlchemy
def db_citations_retrieval(user_input):
	result = engine.execute('select title, author, journal, pubdate, url from citations where citesPmid = :0', [user_input])
	apa_citations = []
	db_journals = []
	db_dates = []
	db_urls = []
	for row in result:
		title = row["title"]
		author = row["author"]

		journal = row["journal"]
		db_journals.append(journal)

		pubdate = row["pubdate"]
		db_dates.append(pubdate)

		url = row["url"]
		db_urls.append(url)

		apa = str(author+' ('+pubdate+'). '+title+'. '+journal+'. Retrieved from '+url)
		apa_citations.append(apa)
	return apa_citations, db_journals, db_dates, db_urls



#Input: pmid that is cited
#Output: list of titles for citation_venn.py
#updated to sqlAlchemy
def db_citation_titles(user_input):
	result = engine.execute('select title from citations where citesPmid = :0', [user_input])
	db_titles = []
	for row in result:
		title = row["title"]
		db_titles.append(title)
	return db_titles


#Input: pmid that is cited
#Output: list of urls for heatmap and barchart hrefs
#Updated to sqlAlchemy
def db_citation_urls(user_input):
	result = engine.execute('select url from citations where citesPmid = :0', [user_input])
	db_urls = []
	for row in result:
		url = row["url"]
		db_urls.append(url)
	return db_urls


#Input: pmid that is cited
#Output: list of pmc_ids for citation
#Updated to sqlalchemy
def db_citation_pmc_ids(user_input):
	result = engine.execute('select pmcid from citations where citesPmid = :0', [user_input])
	db_pmcids = []
	for row in result:
		pmcid = row["pmcid"]
		db_pmcids.append(pmcid)
	# TODO: May 23, 2017 this function PREVIOUSLY did NOT work with pmcidAnnotated() function IF the connection + cursor are closed here
	# TODO: but there were some locked db problems.... Not sure if fixed with sqlAlchemy accessing
	return db_pmcids


#Input: pmid that is cited
#Output: journals and dates of all citing pmcids
#Updated to sqlAlchemy
def db_journals_and_dates(pmid):
	journals = []
	dates = []
	pmcids = []
	result = engine.execute('select pmcid, journal, pubdate from citations where citesPmid = :0', [pmid])
	for row in result:
		pmc = row["pmcid"]
		pmcids.append(pmc)
		j = row["journal"]
		journals.append(j)
		d = row["pubdate"]
		dates.append(d)
	return pmcids, journals, dates



#this function will tell if a pmcid already exists in 'citations' so that it
#doesn't need to be downloaded, scraped, & annotated in Entrez_IR.py and multipreprocessing.py
#input: pmcid
#output: if the pmcid exists in the db, get the entry to copy for the new citing document
#Updated to sqlAlchemy
def checkForPMCID(citation):
	try:
		result = engine.execute('select * from citations where pmcid = :0', [citation])
		exist = result.fetchall()
		#if the row does NOT exist
		if len(exist) == 0:
			record = 'empty'
		#if the row does exist
		else:
			record = exist
	except Exception as e:
		record = 'empty'
	return record


#This function will check if a pmcid has been scraped before, by seeing if anything is in abstract, whole_article, sents, tokens
#Input pmcid, pmid
#Output: if abstract, whole_article etc have fields, return those
#Output: else, result should return as a string "empty"
#Updated to sqlalchemy
def checkIfScraped(citation, user_input):
	s = select([citations.c.abstract, citations.c.whole_article]).\
		where(citations.c.pmcid == citation).\
		where(citations.c.citesPmid == user_input).\
		where(citations.c.abstract != None)
	conn = connection()
	c = conn.execute(s)
	exist = c.fetchone()
	if exist is None: #if the abstract, whole_article columns are blank, its empty
		result = 'empty'
	else:
		result = 'occupied' #otherwise its been scraped
	return result



#Check if a pmcid has been annotated before
#Input: a pmcid (string)
#Output: if the pmcid is in the db AND annotated --> 'yes'
#Output: if the pmcid is in the db AND not annotated --> 'no'
#Output: if the pmcid is NOT in the db at all --> 'empty'
#Updated to sqlAlchemy
def pmcidAnnotated(pmcid):
	try:
		result = engine.execute('select annotated from citations where pmcid = :0', [pmcid])
		#c.execute('''SELECT annotated FROM citations WHERE pmcid=?''', (pmcid,))
		exist = result.fetchone()
		# if the row does NOT exist
		if exist is None:
			record = 'empty'
		# if the row does exist
		else:
			answer = exist[0]
			if answer == 'yes':
				record = 'yes'
			elif answer == 'no':
				record = 'no'
			else:
				record = 'empty'
	except Exception as e:
		record = 'empty'
	return record


#Updated to sqlalchemy
def update_annotations(b, user_input):
	pc = str(b["pmcid"])
	s = b["num_sentences"][0]
	t = b["num_tokens"][0]
	up = citations.update().\
		where(citations.c.pmcid == pc).\
		where(citations.c.citesPmid == user_input).\
		values(dict(sents = s, tokens = t))
	conn = connection()
	conn.execute(up)


######################## SUPPORT FUNCTIONS FOR queries TABLE #######################################
#check if a query is in db
#Updated to sqlAlchemy
def checkForQuery(query):
	result = engine.execute('select query from queries where query= :0', [query])
	exist = result.fetchone()
	if exist is None:
		record = 'empty'
	else:
		record = 'yes'
	return record


#Get information for cached journal data vis
#Updated to sqlalchemy
def getJournalsVis(query):
	result = engine.execute('select range_years, unique_pubs, unique_journals from queries where query = :0', [query])
	for row in result:
		range_years = row["range_years"]
		unique_pubs = row["unique_pubs"]
		unique_journals = row["unique_journals"]
	return range_years, unique_pubs, unique_journals



#################### SUPPORT FUNCTIONS FOR annotations TABLE ############

#checks whether or not a pmcid is in the db
#updated to sqlalchemy
def annotationsCheckPmcid(pmcid):
	result = engine.execute('select pmcid from annotations where pmcid= :0', [pmcid])
	exist = result.fetchone()
	if exist is None:
		record = 'empty'
	else:
		record = 'yes'
	return record


#retrieve the lemmas for citing documents as list of strings
#data should be a list of strings for the documents
#Updated to sqlAlchemy
def getDataSamples(pmcid_list):
	data_samples = []
	pmcid_set = set(pmcid_list) #we only want UNIQUE pmcids
	for pmcid in pmcid_set:
		result = engine.execute('select lemmas from annotations where pmcid = :0', [pmcid])
		#c.execute('''SELECT lemmas FROM annotations WHERE pmcid=?''', (pmcid,))
		for row in result:
			lemmas = row["lemmas"]
			data_samples.append(lemmas)
	return data_samples, pmcid_set




