from flask import Flask
import time, datetime, re
from collections import defaultdict
from sqlalchemy import create_engine, MetaData, Table, select
import logging


'''
Tables in pmids_info.db
1) inputPapers: stores info about the input pmids and their citations
2) citations: stores information about pmcids that cite input pmids
3) queries: stores information about a particular input query e.g. paper1+paper2
4) annotations: stores information about the annotated pmcids' lemmas and named entities categories

'''

def connect_db():
	app = Flask(__name__, static_url_path='/hclent/Webdev-for-bioNLP-lit-tool/flask/static')
	app.config.from_pyfile('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/configscke.cfg', silent=False)  # pass abs path
	engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'])
	return engine


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

logging.basicConfig(filename='.app.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

engine = connect_db()
conn = connection()
inputPapers, citations, queries, annotations = load_tables()


#################### SUPPORT FUNCTIONS FOR inputPapers TABLE ######################
#Input: pmid
#Output: apa citation of *THAT* pmid
#Updated to use sqlalchemy
def db_inputPapers_retrieval(user_input):
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
	return apa


#Input: pmid
#Output: number of citations
#Updated to use sqlalchemy
def db_input_citations_count(user_input):
	s = select([inputPapers.c.num_citations]).\
		where(inputPapers.c.pmid == user_input)
	result = conn.execute(s)
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
	return pmidDict, pmcDict


# After you've scraped the input_paper, write the information about abstract, article, and self_pmcid to db
# from function getSelfText(user_input)
# puts pmcid, abstract check, and article check into db from
# Updated to sqlalchemy syntax
def updateInputPapers(user_input, self_pmcid, abstract, article):
	up = inputPapers.update().\
		where(inputPapers.c.pmid == user_input).\
		values(dict(abstract=abstract, whole_article=article, pmcid=self_pmcid))
	conn.execute(up)



# convert pmid2pmcid
# updated to sqlAlchemy
def pmid2pmcid(user_input):
	s = select([inputPapers.c.pmcid]).\
		where(inputPapers.c.pmid == user_input)
	result = conn.execute(s)
	for row in result:
		pmcid = row["pmcid"]  # return first thing in tuple ('2836516',)
	return pmcid
	# will return NoneType if its not there apparently :)


#This is for TextCompare x-axis labels for eligible_papers (inputPapers)
def db_pmid_axis_label(pmid):
	s = select([inputPapers.c.author, inputPapers.c.pubdate, inputPapers.c.url]).\
		where(inputPapers.c.pmid == pmid)
	result = conn.execute(s)
	label = []
	href_label = []
	for row in result:
		author = row["author"] #grab first author's last name (ie texts up until first space)
		first_author = re.match('^\w+\s' ,author)
		if first_author:
			keep_author = str(first_author.group(0))
		else:
			keep_author = str(author[:5]) + '...'

		pubdate = row["pubdate"] #need to just get year
		year = re.search('\d{4}', pubdate)
		if year:
			keep_year = int(year.group(0))
		else:
			try:
				keep_year = str(pubdate)
			except Exception as e:
				#we need an int no matter what...
				keep_year = "0000"

		url = row["url"]
		apa = str(str(keep_year) + ', ' + keep_author) #2008, author last name
		label.append(apa)
	return label


def db_pmid_hyperlink_retrieval(pmid):
	s = select([inputPapers.c.title, inputPapers.c.author, inputPapers.c.journal, inputPapers.c.pubdate, inputPapers.c.url]).\
		where(inputPapers.c.pmid == pmid)
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
	return apa_citations

######################### SUPPORT FUNCTIONS FOR citations TABLE ###########################
#Input: pmcid
#Output: apa citations for citing pmCids as hyperlinks
#Updated to sqlAlchemy
#UPDATED
def db_citations_hyperlink_retrieval(pmcid):
	s = select([citations.c.title, citations.c.author, citations.c.journal, citations.c.pubdate, citations.c.url]).\
		where(citations.c.pmcid == pmcid)
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
	return apa_citations




#Gets the YEAR only for a given pmcid for the heatmap vis
def db_citations_mini_year(pmcid):
	s = select([citations.c.pubdate]).\
		where(citations.c.pmcid == pmcid)
	result = conn.execute(s)
	years = []
	for row in result:
		pubdate = row["pubdate"] #need to just get year
		year = re.search('\d{4}', pubdate)
		if year:
			keep_year = int(year.group(0))
		else:
			try:
				keep_year = str(pubdate)
			except Exception as e:
				#we need an int no matter what...
				keep_year = "0000"
		years.append(keep_year) #there could be more than one result
	use_year = years[0] #but we'll just use the first result
	return use_year


#This is for the heatmap x-axis labels and TextCompare x-axis labels
def db_citations_mini_hyperlink(pmcid):
	s = select([citations.c.author, citations.c.pubdate, citations.c.url]).\
		where(citations.c.pmcid == pmcid)
	result = conn.execute(s)
	label = []
	for row in result:
		author = row["author"] #grab first author's last name (ie texts up until first space)
		first_author = re.match('^\w+\s' ,author)
		if first_author:
			keep_author = str(first_author.group(0))
		else:
			keep_author = str(author[:5]) + '...'

		pubdate = row["pubdate"] #need to just get year
		year = re.search('\d{4}', pubdate)
		if year:
			keep_year = int(year.group(0))
		else:
			try:
				keep_year = str(pubdate)
			except Exception as e:
				#we need an int no matter what...
				keep_year = "0000"

		url = row["url"]
		apa = str(str(keep_year) + ', ' + keep_author) #2008, author last name
		#Use the href_label if they should be hyperlinks
		#href_label = str('<a href="'+url+'">'+str(apa)+'</a>')
		label.append(apa)
	return label



#Input: pmid
#Output: list of apa citations of pmc-ids citing that pmid
#Updated to sqlAlchemy
def db_citations_retrieval(user_input):
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
	return apa_citations, db_journals, db_dates, db_urls


#Input: pmid that is cited
#Output: list of titles for citation_venn.py
#updated to sqlAlchemy
def db_citation_titles(user_input):
	s = select([citations.c.title]).\
		where(citations.c.citesPmid == user_input)
	result = conn.execute(s)
	db_titles = []
	for row in result:
		title = row["title"]
		db_titles.append(title)
	return db_titles


#Input: pmid that is cited
#Output: list of urls for heatmap and barchart hrefs
#Updated to sqlAlchemy
def db_citation_urls(user_input):
	s = select([citations.c.url]).\
		where(citations.c.citesPmid == user_input)
	result = conn.execute(s)
	db_urls = []
	for row in result:
		url = row["url"]
		db_urls.append(url)
	return db_urls


#Input: pmid that is cited
#Output: list of pmc_ids for citation
#Updated to sqlalchemy
def db_citation_pmc_ids(user_input):
	s = select([citations.c.pmcid]).\
		where(citations.c.citesPmid == user_input)
	result = conn.execute(s)
	db_pmcids = []
	for row in result:
		pmcid = row["pmcid"]
		db_pmcids.append(pmcid)
	return db_pmcids


#Input: pmid that is cited
#Output: journals and dates of all citing pmcids
#Updated to sqlAlchemy
def db_journals_and_dates(pmid):
	journals = []
	dates = []
	pmcids = []
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
	return pmcids, journals, dates



#this function will tell if a pmcid already exists in 'citations' so that it
#doesn't need to be downloaded, scraped, & annotated in Entrez_IR.py and multipreprocessing.py
#input: pmcid
#output: if the pmcid exists in the db, get the entry to copy for the new citing document
#Updated to sqlAlchemy
def checkForPMCID(citation):
	try:
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
	except Exception as e:
		record = 'empty'
	return record


#Updated to sqlalchemy
def update_annotations(b, user_input):
	pc = str(b["pmcid"])
	s = b["num_sentences"]
	t = b["num_tokens"][0]
	up = citations.update().\
		where(citations.c.pmcid == pc).\
		where(citations.c.citesPmid == user_input).\
		values(dict(sents = s, tokens = t))
	conn.execute(up)


######################## SUPPORT FUNCTIONS FOR queries TABLE #######################################
#check if a query is in db
#Updated to sqlAlchemy
def checkForQuery(query):
	s = select([queries.c.query]).\
		where(queries.c.query == query)
	result = conn.execute(s)
	exist = result.fetchone()
	if exist is None:
		record = 'empty'
	else:
		record = 'yes'
	return record


#Get information for cached journal data vis
#Updated to sqlalchemy
def getJournalsVis(query):
	s = select([queries.c.range_years, queries.c.unique_pubs, queries.c.unique_journals]).\
		where(queries.c.query == query)
	result = conn.execute(s)
	for row in result:
		range_years = row["range_years"]
		unique_pubs = row["unique_pubs"]
		unique_journals = row["unique_journals"]
	return range_years, unique_pubs, unique_journals


#################### SUPPORT FUNCTIONS FOR annotations TABLE ############

#checks whether or not a pmcid is in the db
#updated to sqlalchemy
def annotationsCheckPmcid(pmcid):
	s = select([annotations.c.pmcid]).\
		where(annotations.c.pmcid == pmcid)
	result = conn.execute(s)
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
		s = select([annotations.c.lemmas]).\
			where(annotations.c.pmcid == pmcid)
		result = conn.execute(s)
		for row in result:
			lemmas = row["lemmas"]
			data_samples.append(lemmas)
	return data_samples, pmcid_set

########### O R M ########################




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

