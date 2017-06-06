import time, datetime
from collections import defaultdict
from sqlalchemy import create_engine, MetaData, Table, select

#May 26, 2017: accessing db through sqlite3 library depreciated

'''
Tables in pmids_info.db
1) inputPapers: stores info about the input pmids and their citations
2) citations: stores information about pmcids that cite input pmids
3) queries: stores information about a particular input query e.g. paper1+paper2
4) annotations: stores information about the annotated pmcids' lemmas and named entities categories

'''

#Connect to db through SqlAlchemy:
#engine = create_engine('sqlite:///pmids_info.db')
#TODO: DO NOT COMMIT DB PASSWORD TO REPO!!!!! :P
engine = create_engine("mysql://info/for/db?charset=utf8mb4")

def connection():
	conn = engine.connect()
	return conn


#TODO: should conn be a global variable??

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
	conn = connection()
	s = select([inputPapers.c.title, inputPapers.c.author, inputPapers.c.journal, inputPapers.c.pubdate, inputPapers.c.url]).\
		where(inputPapers.c.pmid == user_input)
	result = conn.execute(s)
	for row in result:
		title = row['title']
		author = row['author']
		journal = row['journal']
		pubdate = row['pubdate']
		url = row['url']
		apa = str(author+' ('+pubdate+'). '+title+'. '+journal+'. Retrieved from '+url)
	result.close()
	return apa


#Input: pmid
#Output: number of citations
#Updated to use sqlalchemy
def db_input_citations_count(user_input):
	conn = connection()
	s = select([inputPapers.c.num_citations]).\
		where(inputPapers.c.pmid == user_input)
	result = conn.execute(s)
	for row in result:
		num_citations = row['num_citations']
	result.close()
	return num_citations


# Input: pmid that is cited
# Output: dicts needed for statistics tab
# Updated to use sqlalchemy
def db_statistics(user_input):
	pmidDict = defaultdict(int)
	pmcDict = defaultdict(list)

	num_citations = db_input_citations_count(user_input)
	pmidDict[user_input] = num_citations

	conn = connection()
	s = select([citations.c.pmcid, citations.c.abstract, citations.c.whole_article, citations.c.sents, citations.c.tokens]).\
		where(citations.c.citesPmid == user_input)
	result = conn.execute(s)

	for row in result:
		pmcid = row['pmcid']
		abstract = row['abstract']
		whole = row['whole_article']
		sents = row['sents']
		tokens = row['tokens']
		pmcDict[pmcid] = [abstract, whole, sents, tokens]
	result.close()
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
	conn = connection()
	s = select([inputPapers.c.pmcid]).\
		where(inputPapers.c.pmid == user_input)
	result = conn.execute(s)
	for row in result:
		pmcid = row["pmcid"]  # return first thing in tuple ('2836516',)
	result.close()
	return pmcid
	# will return NoneType if its not there apparently :)



######################### SUPPORT FUNCTIONS FOR citations TABLE ###########################
#Input: pmid
#Output: apa citations for citing pmCids as hyperlinks
#Updated to sqlAlchemy
def db_citations_hyperlink_retrieval(user_input):
	conn = connection()
	s = select([citations.c.title, citations.c.author, citations.c.journal, citations.c.pubdate, citations.c.url]).\
		where(citations.c.citesPmid == user_input)
	result = conn.execute(s)
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
	result.close()
	return apa_citations


#Input: pmid
#Output: list of apa citations of pmc-ids citing that pmid
#Updated to sqlAlchemy
def db_citations_retrieval(user_input):
	conn = connection()
	s = select([citations.c.title, citations.c.author, citations.c.journal, citations.c.pubdate, citations.c.url]).\
		where(citations.c.citesPmid == user_input)
	result = conn.execute(s)
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
	result.close()
	return apa_citations, db_journals, db_dates, db_urls


#Input: pmid that is cited
#Output: list of titles for citation_venn.py
#updated to sqlAlchemy
def db_citation_titles(user_input):
	conn = connection()
	s = select([citations.c.title]).\
		where(citations.c.citesPmid == user_input)
	result = conn.execute(s)
	db_titles = []
	for row in result:
		title = row["title"]
		db_titles.append(title)
	result.close()
	return db_titles


#Input: pmid that is cited
#Output: list of urls for heatmap and barchart hrefs
#Updated to sqlAlchemy
def db_citation_urls(user_input):
	conn = connection()
	s = select([citations.c.url]).\
		where(citations.c.citesPmid == user_input)
	result = conn.execute(s)
	db_urls = []
	for row in result:
		url = row["url"]
		db_urls.append(url)
	result.close()
	return db_urls


#Input: pmid that is cited
#Output: list of pmc_ids for citation
#Updated to sqlalchemy
def db_citation_pmc_ids(user_input):
	conn = connection()
	s = select([citations.c.pmcid]).\
		where(citations.c.citesPmid == user_input)
	result = conn.execute(s)
	db_pmcids = []
	for row in result:
		pmcid = row["pmcid"]
		db_pmcids.append(pmcid)
	result.close()
	return db_pmcids


#Input: pmid that is cited
#Output: journals and dates of all citing pmcids
#Updated to sqlAlchemy
def db_journals_and_dates(pmid):
	journals = []
	dates = []
	pmcids = []
	conn = connection()
	s = select([citations.c.pmcid, citations.c.journal, citations.c.pubdate]).\
		where(citations.c.citesPmid == pmid)
	result = conn.execute(s)
	for row in result:
		pmc = row["pmcid"]
		pmcids.append(pmc)
		j = row["journal"]
		journals.append(j)
		d = row["pubdate"]
		dates.append(d)
	result.close()
	return pmcids, journals, dates



#this function will tell if a pmcid already exists in 'citations' so that it
#doesn't need to be downloaded, scraped, & annotated in Entrez_IR.py and multipreprocessing.py
#input: pmcid
#output: if the pmcid exists in the db, get the entry to copy for the new citing document
#Updated to sqlAlchemy
def checkForPMCID(citation):
	try:
		conn = connection()
		s = citations.select().\
			where(citations.c.pmcid == citation)
		result = conn.execute(s)
		exist = result.fetchall()
		#if the row does NOT exist
		if len(exist) == 0:
			record = 'empty'
		#if the row does exist
		else:
			record = exist
		result.close()
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
		conn = connection()
		s = select([citations.c.annotated]).\
			where(citations.c.pmcid == pmcid)
		result = conn.execute(s)
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
		result.close()
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
	conn = connection()
	s = select([queries.c.query]).\
		where(queries.c.query == query)
	result = conn.execute(s)
	exist = result.fetchone()
	if exist is None:
		record = 'empty'
	else:
		record = 'yes'
	result.close()
	return record


#Get information for cached journal data vis
#Updated to sqlalchemy
def getJournalsVis(query):
	conn = connection()
	s = select([queries.c.range_years, queries.c.unique_pubs, queries.c.unique_journals]).\
		where(queries.c.query == query)
	result = conn.execute(s)
	for row in result:
		range_years = row["range_years"]
		unique_pubs = row["unique_pubs"]
		unique_journals = row["unique_journals"]
	result.close()
	return range_years, unique_pubs, unique_journals



#################### SUPPORT FUNCTIONS FOR annotations TABLE ############

#checks whether or not a pmcid is in the db
#updated to sqlalchemy
def annotationsCheckPmcid(pmcid):
	conn = connection()
	s = select([annotations.c.pmcid]).\
		where(annotations.c.pmcid == pmcid)
	result = conn.execute(s)
	exist = result.fetchone()
	if exist is None:
		record = 'empty'
	else:
		record = 'yes'
	result.close()
	return record


#retrieve the lemmas for citing documents as list of strings
#data should be a list of strings for the documents
#Updated to sqlAlchemy
def getDataSamples(pmcid_list):
	data_samples = []
	pmcid_set = set(pmcid_list) #we only want UNIQUE pmcids
	conn = connection()
	for pmcid in pmcid_set:
		s = select([annotations.c.lemmas]).\
			where(annotations.c.pmcid == pmcid)
		result = conn.execute(s)
		for row in result:
			lemmas = row["lemmas"]
			data_samples.append(lemmas)
	result.close()
	return data_samples, pmcid_set


'''
Database setup for MySql:
NOTE: I did not have the foresight to initialize the tables with utf8mb4 encoding,
      thus I had to go back and
      ALTER DATABASE scke CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
      & for each table
       ALTER TABLE <tablename> CONVERT TO CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;

CREATE TABLE inputPapers (post_id MEDIUMINT NOT NULL AUTO_INCREMENT, datestamp DATE, pmid VARCHAR(15), title VARCHAR(1000), author VARCHAR(10000), journal VARCHAR(500), pubdate VARCHAR(50), url VARCHAR(200), num_citations INT, abstract VARCHAR(10), whole_article VARCHAR(10), pmcid VARCHAR(20), KEY(post_id));

CREATE TABLE citations (post_id MEDIUMINT NOT NULL AUTO_INCREMENT, datestamp DATE, pmcid VARCHAR(15), title VARCHAR(1000), author VARCHAR(10000), journal VARCHAR(500), pubdate VARCHAR(50), citesPmid VARCHAR(20), url VARCHAR(200), abstract VARCHAR(10), whole_article VARCHAR(10), sents INT, tokens INT, annotated VARCHAR(10), KEY(post_id));

CREATE TABLE queries (post_id MEDIUMINT NOT NULL AUTO_INCREMENT, datestamp DATE, query VARCHAR(500), range_years VARCHAR(20), total_pubs INT, unique_pubs INT, unique_journals INT, num_abstracts INT, num_whole_articles INT, num_sents INT, num_tokens INT, KEY(post_id));

CREATE TABLE annotations (post_id MEDIUMINT NOT NULL AUTO_INCREMENT, datestamp DATE, pmcid VARCHAR(20), lemmas LONGTEXT, bioprocess LONGTEXT, cell_lines LONGTEXT, cell_components LONGTEXT, family LONGTEXT, gene_product LONGTEXT, organ LONGTEXT, simple_chemical LONGTEXT, site LONGTEXT, species LONGTEXT, tissue_type LONGTEXT, KEY(post_id));

'''

