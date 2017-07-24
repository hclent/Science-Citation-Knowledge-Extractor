from __future__ import print_function
import time, sys, pickle, datetime, os.path, logging, json, re
import xml.etree.ElementTree as ET
from time import sleep
from Bio import Entrez
from configapp import app
from database_management import checkForPMCID, checkIfScraped


#Entrez Information Retrieval
#This code uses BioPython to access NCBI's API (Entrez)
#From the NCBI API, this code references PubMed and PubMedCentral
#With an input PubMed ID's (pmids), the code will retrieve information about the original pmid,
#and information about pubplications that cite this pmid via PubMedCentral ID's (pmcids)


Entrez.email = "hclent1@gmail.com" 
Entrez.tool = "MyInfoRetrieval"

#Create log
logging.basicConfig(filename='.app.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logging.info('Started')


#Input: pmid
#Output: self info
#[('KIFC1-like motor protein associates with the cephalopod manchette and participates in sperm nuclear morphogenesis in Octopus tankahkeei.',
# ['Wang W', 'Zhu JQ', 'Yu HM', 'Tan FQ', 'Yang WX'],
# 'PloS one',
# '2010 Dec 20',
# 'www.blah')]
def getMainInfo(pmid):
	logging.info("beginning function getMainInfo")
	t0 = time.time()
	logging.info("making handle ... ")
	handle = Entrez.esummary(db="pubmed", id=pmid)
	logging.info("made the handle!")
	record = Entrez.read(handle)
	logging.info("made a record")
	title = [record[0]["Title"]] #make a list
	logging.info(title)
	authors = [record[0]["AuthorList"]]
	logging.info(authors)
	journal = [record[0]["FullJournalName"]]
	logging.info(journal)
	pubdate = [record[0]["PubDate"]]
	logging.info(pubdate)
	url = ["https://www.ncbi.nlm.nih.gov/pubmed/"+pmid]
	logging.info(url)
	self_info = list(zip(title, authors, journal, pubdate, url))
	logging.info(self_info)
	logging.info("self info: done in %0.3fs." % (time.time() - t0))
	return self_info


#gives pmcid for the pmid
def getAlternativeId(pmid):
	logging.info("in function get AlternativeID")
	try:
		handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
		logging.info("made the handle")
		record = Entrez.read(handle)
		logging.info("made the record")
		data = record[0]["PubmedData"]
		ids = data["ArticleIdList"]
		logging.info(ids)
		pmc_id = ids[3] #PMC12345
		pmcid = re.sub('PMC', '', pmc_id) #12345 without "PMC"
		logging.info("pmcid for input paper" + str(pmcid))
		return pmcid
	except Exception as e:
		logging.info(e)
		logging.info("there was no PMCID for this paper")
		fake_pmcid = "NA"
		return fake_pmcid



#Input: Pmid
#Output: list of pmcids of the articles citing you
def getCitationIDs(pmid): #about the same speed as MainCrawl.py
	t0 = time.time()
	results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", from_uid=pmid))
	pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
	logging.info("get Citation PMCIDs: done in %0.3fs." % (time.time() - t0))
	time.sleep(3) #sleep for 3 seconds to give the API a rest
	return pmc_ids #list


#Input: a single citation (pmcid)
#Output: title, authors, journals, date, url
def connectToNCBI(citation):
	handle = Entrez.esummary(db="pmc", id=citation)
	record = Entrez.read(handle)
	t = record[0]["Title"]
	a = record[0]["AuthorList"]
	j = record[0]["FullJournalName"]
	d = record[0]["PubDate"]
	u = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC"+citation
	return t, a, j, d, u


#Input: Citing pmcids
#Output: Basic info about these pmcids -- list of dictionaries
#Update: Changed from list of lists to dictionary so that the code doesn't rely on indexing
def getCitedInfo(pmcid_list, pmid, conn):
	t0 = time.time()
	i = 1
	allCitations = []
	for citation in pmcid_list:
		logging.info("citation no. " + str(i) + " ...")
		#first check the db to see if this citation already exists
		rows = checkForPMCID(citation, conn) #Now returns all rows with that pmcid (if any)
		if rows != 'empty':
			logging.info("the pmcid already exists in db")
			# if row(s) exists, need to check what pmids they are citing
			cited_pmids_in_db = [r["citesPmid"] for r in rows]
			if pmid in cited_pmids_in_db:
				#if our pmid is in there, we don't need to copy the entry
				logging.info("doesnt need to be coppied again")
			#if our pmid isn't in there, we do need to copy it
			if pmid not in cited_pmids_in_db:
				logging.info("pmcid needs to be coppied")
				row = rows[0] #just grab the first one
				pmcid = citation
				title = row["title"]
				author = row["author"]
				journal = row["journal"]
				date = row["datestamp"]
				citesPmid = pmid #the input pmid
				other_citesPmid = row["citesPmid"] #the pmid from citesPmid may or may not be the same as the input parameter pmid
				url = row["url"]
				abstract_check = row["abstract"] #None or yes/no
				article_check = row["whole_article"] #None or yes/no
				sents = row["sents"] #None or number
				tokens = row["tokens"] #None or number
				annotated = row["annotated"] #None or yes/no

				#Cannot write to DB in this function because the db will be locked, as it is simultaneously reading/writing.
				#return a more-complete citationsDict to content_management to be written to the db there
				citationsDict = {'citesPmid': pmid, "pmcid": citation, "pmc_titles": [], "pmc_authors": [],
							 "pmc_journals": [], "pmc_dates": [], "pmc_urls": [], "abstract_check": [], "article_check":[],
							 "sents": [], "tokens":[], "annotated":[]}
				citationsDict["pmc_titles"].append(title)
				citationsDict["pmc_authors"].append(author)
				citationsDict["pmc_journals"].append(journal)
				citationsDict["pmc_dates"].append(date)
				citationsDict["pmc_urls"].append(url)
				citationsDict["abstract_check"].append(abstract_check)
				citationsDict["article_check"].append(article_check)
				citationsDict["sents"].append(sents)
				citationsDict["tokens"].append(tokens)
				citationsDict["annotated"].append(annotated)
				allCitations.append(citationsDict)

		#else, if the pmcid citation has never been seen before by the db, retrieve info
		else:
			logging.info("the pmcid has never been seen by the db before")
			citationsDict = {'citesPmid': pmid, "pmcid": citation, "pmc_titles": [], "pmc_authors": [], "pmc_journals": [],
							 "pmc_dates": [], "pmc_urls": []}
			#sometimes the connection fails: Need to re-try if this happens
			try:
				t, a, j, d, u = connectToNCBI(citation)
				citationsDict["pmc_titles"].append(t)
				citationsDict["pmc_authors"].append(a)
				citationsDict["pmc_journals"].append(j)
				citationsDict["pmc_dates"].append(d)
				citationsDict["pmc_urls"].append(u)
			except Exception as e:
				logging.info(str(e))
				logging.info("Failed to connect to NCBI. Lets try again in 3 seconds... Retry 1/3")
				time.sleep(3)
				t, a, j, d, u = connectToNCBI(citation)
				citationsDict["pmc_titles"].append(t)
				citationsDict["pmc_authors"].append(a)
				citationsDict["pmc_journals"].append(j)
				citationsDict["pmc_dates"].append(d)
				citationsDict["pmc_urls"].append(u)
			except Exception as e2:
				logging.info(str(e2))
				logging.info("Failed to connect to NCBI again. Let's try again in 5 seconds ... Retry 2/3")
				time.sleep(5)
				t, a, j, d, u = connectToNCBI(citation)
				citationsDict["pmc_titles"].append(t)
				citationsDict["pmc_authors"].append(a)
				citationsDict["pmc_journals"].append(j)
				citationsDict["pmc_dates"].append(d)
				citationsDict["pmc_urls"].append(u)
			except Exception as e3:
				logging.info(str(e3))
				logging.info("Failed to connect to NCBI a third time. Try once more in 8 seconds before giving up ... Retry 3/3")
				time.sleep(5)
				t, a, j, d, u = connectToNCBI(citation)
				citationsDict["pmc_titles"].append(t)
				citationsDict["pmc_authors"].append(a)
				citationsDict["pmc_journals"].append(j)
				citationsDict["pmc_dates"].append(d)
				citationsDict["pmc_urls"].append(u)
			except Exception as e4:
				logging.info(str(e4))
				logging.info("Failed to connect to NCBI 4 times. Let's skip this entry :( ")
				pass
			allCitations.append(citationsDict)
			time.sleep(3)
		i += 1
	logging.info("get citations info: done in %0.3fs." % (time.time() - t0))
	return allCitations



#Input: XML string of PMC entry generated with getContentPMC
#Output: Abstract and journal text
def parsePMC(xml_string):
	main_text = []
	root = ET.fromstring(xml_string) #parsed
	abstract_check = [] #did it get the abstract? yes/no?
	whole_article_check = [] #did it get the whole text? yes/no?

	#Get abstract and add to doc
	try:
		#check for new line!
		abstract = root.find('.//abstract')
		full_abs = ("".join(abstract.itertext()))
		logging.info("* Got abstract")
		main_text.append(full_abs)
		if full_abs: #if the abstract exists
			abstract_check.append('yes')
	except Exception as e:
		logging.info("The following PMCID has no abstract")
		string = "Some data"
		main_text.append(string)
		abstract_check.append('no')
	try:
		#Get main text and add to doc
		text = root.findall('.//p')
		#problem: will return an empty list if there is no text :(
		#error handling here:
		if not text:
			logging.info("main text is empty!!!!!!!!!")
			whole_article_check.append('no')
		#if the text really does exist:
		if text:
			whole_article_check.append('yes')
			for t in text:
				full_text = ("".join(t.itertext()))
				main_text.append(full_text)
			logging.info("* Got main text")
	except Exception as e:
		logging.info("The following PMCID has no main text")
		string = "data"
		main_text.append(string)
		whole_article_check.append('no')
	logging.info("ABSTRACT CHECK")
	logging.info(abstract_check)
	logging.info("WHOLE_ARTICLE CHECK")
	logging.info(whole_article_check)
	return main_text, abstract_check, whole_article_check



#Input: the list of pmcids citing some pmid
#For each citing pmc_id, this function gest the xml, which is then parsed by parsePMC()
#Output: Those texts for each pmcid saved as pmcid.txt to the folder pmcid[:3]/pmicd[3:6] for better organization
def getContentPMC(pmcids_list, pmid, conn):
	t0 = time.time()
	i = 1

	contentDictList = []

	for citation in pmcids_list:
		#print(citation)
		# if the pmc is already in the database for another pmid, then don't rescrape, but DO add to
		# a record that this pmc is cited by the new input paper (will be done at this point though... I think?)
		row = checkIfScraped(citation, pmid, conn)
		if row == 'empty':

			# if the pmcid in the database has our query pmid AND has no abstract check and article check,
			# get the XML record, parse the xml, and do the abstract and main text checks
			contentDict = {"pmcid": citation, "citesPmid": pmid, "all_abstract_check": [], "all_article_check": []}
			logging.info("pmcid never seen before. needs to be scraped + annotated ")

			prefix = str(citation[0:3]) #folder for first 3 digits of pmcid
			suffix = str(citation[3:6]) #folder for second 3 digits of pmcid nested in prefix

			try:
				os.makedirs(os.path.join((app.config['PATH_TO_CACHE']), prefix)) # creates folder named after first 3 digits of pmcid
			except OSError:
				if os.path.isdir(os.path.join((app.config['PATH_TO_CACHE']), prefix)):
					pass
				else:
					raise

			try:
				os.makedirs(os.path.join((app.config['PATH_TO_CACHE']), prefix, suffix)) # creates folder named after second 3 digits of pmicd
			except OSError:
				if os.path.isdir(os.path.join((app.config['PATH_TO_CACHE']), prefix, suffix)):
					pass
				else:
					raise

			logging.info(str(i)+" paper")
			logging.info("CITATION: " +str(citation))
			handle = Entrez.efetch(db="pmc", id=citation, rettype='full', retmode="xml")
			xml_record = handle.read() #xml str
			#print(xml_record)
			logging.info("* got xml record")
			main_text, abstract_check, whole_article_check = parsePMC(xml_record)
			for yn in abstract_check: #'yn' = 'yes no'
				contentDict["all_abstract_check"].append(yn)
			for yn in whole_article_check:
				contentDict["all_article_check"].append(yn)
			#append yes/no content dicts to list
			contentDictList.append(contentDict)
			#print
			logging.info("* ready to print it")
			completeName = os.path.join((app.config['PATH_TO_CACHE']), prefix, suffix, (str(citation)+'.txt'))  #pmcid.txt #save to suffix path
			#print(completeName)
			sys.stdout = open(completeName, "w")
			print(main_text)
			i += 1
			time.sleep(3)
		else:
			logging.info("the document isn't new. pass")
			#TODO: check that we really have the document
			pass
	logging.info("got documents: done in %0.3fs." % (time.time() - t0))
	return contentDictList



#retrieve the txt for input papers as well
def getSelfText(pmid):
	pmcid_list = [getAlternativeId(pmid)] #pmid --> list of 1 pmcid
	fake_list = ["NA"] #not avaliable
	if pmcid_list == fake_list:
		logging.info("there wasn't a matching PMCID for the input PMID")
		fake_pmcid = fake_list[0] #str
		fake_abstract = "NA" #str
		fake_article = "NA" #str
		return fake_pmcid, fake_abstract, fake_article #if there's nothing, we'll fill the db with "NA"
	else:
		logging.info(pmcid_list)
		self_pmcid = pmcid_list[0] #str
		self_abstract_check, self_article_check = getContentPMC(pmcid_list)
		abstract = self_abstract_check[0] #str
		article = self_article_check[0] #str
		return self_pmcid, abstract, article


################### Notes ##############
#Rarely, the XML will return this:
	# <?xml version="1.0"?>
	# <!DOCTYPE pmc-articleset PUBLIC "-//NLM//DTD ARTICLE SET 2.0//EN" "http://dtd.nlm.nih.gov/ncbi/pmc/articleset/nlm-articleset-2.0.dtd">
	# <pmc-articleset><Reply Id="24948109" error="The following PMCID is not available: 24948109"/></pmc-articleset>
#Thus pasrePMC() has exception handeling for this

# Beginning August, NCBI is having connection errors: "pmc is not db". Must add exception handling for it