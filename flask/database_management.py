import sqlite3, time, datetime
from collections import defaultdict


#Basic SQLITE3 structure
#Database in webdev-biotool for managing pmids and scraped webpages
conn = sqlite3.connect(database='pmids_info.db',timeout=5) #connect to database
c = conn.cursor() #cursor


#Define connection and cursor
def connection():
	conn = sqlite3.connect(database='pmids_info.db', timeout=5) #connect to database
	c = conn.cursor() #cursor
	return conn, c


#Input: pmid
#Output: apa citation of *THAT* pmid
def db_inputPapers_retrieval(user_input):
	c.execute('''SELECT title, author, journal, pubdate, url FROM inputPapers WHERE pmid=?''', (user_input,))
	for row in c:
		title = row[0]
		author = row[1]
		journal = row[2]
		pubdate = row[3]
		url = row[4]
		apa = str(author+' ('+pubdate+'). '+title+'. '+journal+'. Retrieved from '+url)
		return apa

#Input: pmid
#Output: number of citations
def db_input_citations_count(user_input):
	c.execute('''SELECT num_citations FROM inputPapers WHERE pmid=?''', (user_input,))
	for row in c:
		num_citations = row[0]
		return num_citations


#Input: pmid
#Output: apa citations for citing pmCids as hyperlinks
def db_citations_hyperlink_retrieval(user_input):
	c.execute('''SELECT title, author, journal, pubdate, url FROM citations WHERE citesPmid=?''', (user_input,))
	apa_citations = []
	for row in c:
		title = row[0]
		author = row[1]
		journal = row[2]
		pubdate = row[3]
		url = row[4]
		apa = str(author+' ('+pubdate+'). '+title+'. '+journal+'. Retrieved from '+url)
		href_label = str('<a href="'+url+'">'+str(apa)+'</a>')
		apa_citations.append(href_label)
	return apa_citations


#Input: pmid
#Output: list of apa citations of pmc-ids citing that pmid
def db_citations_retrieval(user_input):
	c.execute('''SELECT title, author, journal, pubdate, url FROM citations WHERE citesPmid=?''', (user_input,))
	apa_citations = []
	db_journals = []
	db_dates = []
	db_urls = []
	for row in c:
		title = row[0]
		author = row[1]

		journal = row[2]
		db_journals.append(journal)

		pubdate = row[3]
		db_dates.append(pubdate)

		url = row[4]
		db_urls.append(url)

		apa = str(author+' ('+pubdate+'). '+title+'. '+journal+'. Retrieved from '+url)
		apa_citations.append(apa)
	return apa_citations, db_journals, db_dates, db_urls


#Input: pmid that is cited
#Output: list of titles for citation_venn.py
def db_citation_titles(user_input):
	c.execute('''SELECT title, author, journal, pubdate, url FROM citations WHERE citesPmid=?''', (user_input,))
	db_titles = []
	for row in c:
		title = row[0]
		db_titles.append(title)
	return db_titles


#Input: pmid that is cited
#Output: list of urls for heatmap and barchart hrefs
def db_citation_urls(user_input):
	c.execute('''SELECT url FROM citations WHERE citesPmid=?''', (user_input,))
	db_urls = []
	for row in c:
		url = row[0]
		db_urls.append(url)
	return db_urls


#Input: pmid that is cited
#Output: list of pmc_ids for citation
def db_citation_pmc_ids(user_input):
	c.execute('''SELECT pmcid FROM citations WHERE citesPmid=?''', (user_input,))
	db_pmcids = []
	for row in c:
		pmcid = row[0]
		db_pmcids.append(pmcid)
	return db_pmcids



#Input: pmid that is cited
#Output: journals and dates of all citing pmcids
def db_journals_and_dates(pmid):
	journals = []
	dates = []
	pmcids = []
	c.execute('''SELECT pmcid, journal, pubdate FROM citations WHERE citesPmid=?''', (pmid,))
	for row in c:
		pmc = row[0]
		pmcids.append(pmc)
		j = row[1]
		journals.append(j)
		d = row[2]
		dates.append(d)
	return pmcids, journals, dates



#Input: pmid that is cited
#output: dicts needed for statistics tab
def db_statistics(user_input):
	pmidDict = defaultdict(int)
	#{pmid: num citations}
	pmcDict = defaultdict(list)
    #dict value [0] = num abstracts
    #dict value [1] = num whole articles
    #dict value [2] = num sentences
    #dict value [3] = num tokens
	c.execute('''SELECT num_citations FROM inputPapers WHERE pmid=?''', (user_input,))
	for row in c:
		total_citations = row[0]
		pmidDict[user_input] = total_citations

	c.execute('''SELECT pmcid, abstract, whole_article, sents, tokens FROM citations WHERE citesPmid=?''', (user_input,))
	for row in c:
		pmcid = row[0]
		abstract = row[1]
		whole = row[2]
		sents = row[3]
		tokens = row[4]
		pmcDict[pmcid] = [abstract, whole, sents, tokens]

	return pmidDict, pmcDict


#After you've scraped the input_paper, write the information about abstract, article, and self_pmcid to db
def updateInputPapers(user_input, self_pmcid, abstract, article):
	#from function getSelfText(user_input)
    #put pmcid, abstract check, and article check into db from
    conn, c = connection()
    c.execute("UPDATE inputPapers SET abstract=?, whole_article=?, pmcid =? WHERE pmid=?",(abstract, article, self_pmcid, user_input))  # put user pmid into db
    conn.commit()


#convert pmid2pmcid
def pmid2pmcid(user_input):
	c.execute('''SELECT pmcid FROM inputPapers WHERE pmid=?''', (user_input,))
	for pmcid in c:
		return pmcid[0] #return first thing in tuple ('2836516',)
		# will return NoneType if its not there apparently :)


def retrieveAllPmcids():
	c.execute('''SELECT pmcid FROM citations''')
	in_db_pmcids = [pmcid[0] for pmcid in c]
	return in_db_pmcids


def retrieveAllPmids():
	c.execute('''SELECT pmid FROM inputPapers''')
	pmid_list = [pmid[0] for pmid in c]
	return pmid_list



#this function will tell if a pmcid already exists in 'citations' so that it
#doesn't need to be downloaded, scraped, & annotated in Entrez_IR.py and multipreprocessing.py
#input: pmcid
#output: if the pmcid exists in the db, get the entry to copy for the new citing document
def checkForPMCID(citation):
	try:
		c.execute('''SELECT pmcid, title, author, journal, pubdate, citesPmid, url, abstract, whole_article,
 					sents, tokens, annotated FROM citations WHERE pmcid=?''',(citation,))
		exist = c.fetchall()
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
def checkIfScraped(citation, user_input):
	c.execute('''SELECT abstract, whole_article FROM citations WHERE pmcid=? AND citesPmid=? AND abstract IS NOT NULL AND whole_article IS NOT NULL''',
			  (citation, user_input,))
	exist = c.fetchone()
	if exist is None: #if the abstract, whole_article columns are blank, its empty
		result = 'empty'
	else:
		result = 'occupied' #otherwise its been scraped
	return result



#Check if a pmcid has been annotated before
def pmcidAnnotated(pmcid):
	try:
		c.execute('''SELECT annotated FROM citations WHERE pmcid=?''', (pmcid,))
		exist = c.fetchone()
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
		print(record)
	return record


#Create table for inputPapers
# def create_table_input():
# 	c.execute('''CREATE TABLE IF NOT EXISTS inputPapers
# 		(post_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
# 		datestamp TEXT,
# 		pmid TEXT,
# 		title TEXT,
# 		author TEXT,
# 		journal TEXT,
# 		pubdate TEXT,
# 		url TEXT,
# 		num_citations NUMBER)''')
#
#

# #Create table for citations of the citations
# def create_table_citations():
# 	c.execute('''CREATE TABLE IF NOT EXISTS citations
# 		(post_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
# 		datestamp TEXT,
# 		pmcid TEXT,
# 		title TEXT,
# 		author TEXT,
# 		journal TEXT,
# 		pubdate TEXT,
# 		citesPmid TEXT,
# 		url TEXT,
# 		abstract TEXT,
# 		whole_article TEXT,
# 		sents NUMBER,
# 		tokens NUMBER)''')
#

#Create table for queries put in to the tool!
# def create_table_citations():
# 	c.execute('''CREATE TABLE IF NOT EXISTS queries
# 		(post_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
# 		datestamp TEXT,
# 		query TEXT,
# 		range_years TEXT,
# 		total_pubs NUMBER,
# 		unique_pubs NUMBER,
# 		unique_journals NUMBER,
# 		num_abstracts NUMBER,
# 		num_whole_articles NUMBER,
# 		num_sents NUMBER,
# 		num_tokens NUMBER)''')
#
# create_table_citations()
# #test
# def test_data_entry():
# 	c.execute("INSERT INTO queries VALUES(1, '05-12-2017', 'aaa+bbb', '1992-2017', 100, 90, 41, 90, 90, 28000, 1000000)")
# 	conn.commit() #to save it to db
#
# 	c.execute("SELECT * FROM queries")
# 	[print(row) for row in c.fetchall()]
#
# 	# c.close()
# 	# conn.close()
#
#
# #print table
# def print_inputPapers():
# 	c.execute("SELECT * FROM queries")
# 	[print(row) for row in c.fetchall()]

