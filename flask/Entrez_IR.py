from __future__ import print_function
from Bio import Entrez
import time, sys, pickle
from time import sleep
import datetime
import os.path, logging, json, re
import xml.etree.ElementTree as ET
import pprint



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
# 'h')]
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
	url = ["http://www.ncbi.nlm.nih.gov/pubmed/"+pmid]
	logging.info(url)
	self_info = list(zip(title, authors, journal, pubdate, url))
	logging.info(self_info)

	logging.info("self info: done in %0.3fs." % (time.time() - t0))
	return self_info

#gives pmcid for the pmid
def getAlternativeId(pmid):
	logging.info("in function get AlternativeID")
	handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
	logging.info("made the handle")
	record = Entrez.read(handle)
	logging.info("made the record")
	data = record[0]["PubmedData"]
	ids = data["ArticleIdList"]
	logging.info(ids)
	try:
		pmc_id = ids[3] #PMC12345
		pmcid = re.sub('PMC', '', pmc_id) #12345 without "PMC"
		logging.info("pmcid for input paper" + str(pmcid))
		#print(pmcid)
		return pmcid
	except Exception as e:
		logging.info(e)
		logging.info("there was no PMCID for this paper")
		fake_pmcid = "None"
		return fake_pmcid


#Input: Pmid
#Output: list of pmcids of the articles citing you
def getCitationIDs(pmid): #about the same speed as MainCrawl.py
	t0 = time.time()
	results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", from_uid=pmid))
	pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
	return pmc_ids #list
	logging.info("get Citation PMCIDs: done in %0.3fs." % (time.time() - t0))
	time.sleep(3)


#Input: a single citation
#Output: title, authors, journals, date, url
def connectToNCBI(citation):
	handle = Entrez.esummary(db="pmc", id=citation)
	record = Entrez.read(handle)
	t = record[0]["Title"]
	a = record[0]["AuthorList"]
	j = record[0]["FullJournalName"]
	d = record[0]["PubDate"]
	u = "http://www.ncbi.nlm.nih.gov/pmc/articles/PMC"+citation
	return t, a, j, d, u


#Input: Citing pmcids
#Output: Basic info about these pmcids
def getCitedInfo(pmcid_list): 
	t0 = time.time()
	pmc_titles = []
	pmc_authors = []
	pmc_journals = []
	pmc_dates = []
	pmc_urls = []
	i = 1
	for citation in pmcid_list:
		logging.info("citation no. " + str(i) + " ...")
		#sometimes the connection fails: Need to re-try
		#
		#Need to re-try if this happens:
		try:
			t, a, j, d, u = connectToNCBI(citation)
			pmc_titles.append(t)
			pmc_authors.append(a)
			pmc_journals.append(j)
			pmc_dates.append(d)
			pmc_urls.append(u)
		except Exception as e:
			logging.info(str(e))
			logging.info("Failed to connect to NCBI. Lets try again in 3 seconds... Retry 1/3")
			time.sleep(3)
			t, a, j, d, u = connectToNCBI(citation)
			pmc_titles.append(t)
			pmc_authors.append(a)
			pmc_journals.append(j)
			pmc_dates.append(d)
			pmc_urls.append(u)
		except Exception as e2:
			logging.info(str(e2))
			logging.info("Failed to connect to NCBI again. Let's try again in 5 seconds ... Retry 2/3")
			time.sleep(5)
			t, a, j, d, u = connectToNCBI(citation)
			pmc_titles.append(t)
			pmc_authors.append(a)
			pmc_journals.append(j)
			pmc_dates.append(d)
			pmc_urls.append(u)
		except Exception as e3:
			logging.info(str(e3))
			logging.info("Failed to connect to NCBI a third time. Try once more in 8 seconds before giving up ... Retry 3/3")
			time.sleep(8)
			t, a, j, d, u = connectToNCBI(citation)
			pmc_titles.append(t)
			pmc_authors.append(a)
			pmc_journals.append(j)
			pmc_dates.append(d)
			pmc_urls.append(u)
		except Exception as e4:
			logging.info(str(e4))
			logging.info("Failed to connect to NCBI 4 times. Let's skip this entry :( ")
			pass
		time.sleep(3)
		i += 1
	#main_info = (list(zip(pmc_titles, pmc_authors, pmc_journals, pmc_urls)))
	logging.info("get citations info: done in %0.3fs." % (time.time() - t0))
	return pmc_titles, pmc_authors, pmc_journals, pmc_dates, pmc_urls


#Input: XML string of PMC entry generated with getContentPMC
#Output: Abstract and journal text
def parsePMC(xml_string, pmid):
	main_text = []
	root = ET.fromstring(xml_string) #parsed
	abstract_check = [] #did it get the abstract? yes/no?
	whole_article_check = [] #did it get the whole text? yes/no?

	#Get abstract and add to doc
	try:
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




#Input: pmid and the list of pmcids citing it
#For each citing pmc_id, this function gest the xml, which is then parsed by parsePMC()
#Output: Journal texts for each pmcid
def getContentPMC(pmid, pmcids_list):
	t0 = time.time()
	i = 1
	dirname = '/home/hclent/data/'+pmid
	all_abstract_check = []
	all_article_check = []
	try:
		os.makedirs(dirname) #creates folder named after pmid
	except OSError:
		if os.path.isdir(dirname):
			pass
		else:
			raise
	for citation in pmcids_list:
		logging.info(str(i)+" paper")
		logging.info("CITATION: " +str(citation))
		handle = Entrez.efetch(db="pmc", id=citation, rettype='full', retmode="xml")
		xml_record = handle.read() #xml str
		#print(xml_record)
		logging.info("* got xml record")
		main_text, abstract_check, whole_article_check = parsePMC(xml_record, pmid)
		for yn in abstract_check:
			all_abstract_check.append(yn)
		for yn in whole_article_check:
			all_article_check.append(yn)
		logging.info("* ready to print it")
		save_path = '/home/hclent/data/'+(str(pmid))+'/' #must save to data, in proper file
		completeName = os.path.join(save_path, (str(pmid)+'_'+str(i)+'.txt'))
		sys.stdout = open(completeName, "w")
		print(main_text)
		i += 1
		time.sleep(3)
	logging.info("got documents: done in %0.3fs." % (time.time() - t0))
	return all_abstract_check, all_article_check


#save self document with prefix "self_12345.txt"
def get_self_ContentPMC(pmid, pmcids_list):
	t0 = time.time()
	i = 1
	dirname = '/home/hclent/data/'+pmid
	all_abstract_check = []
	all_article_check = []
	try:
		os.makedirs(dirname) #creates folder named after pmid
	except OSError:
		if os.path.isdir(dirname):
			pass
		else:
			raise
	for citation in pmcids_list:
		logging.info(str(i)+" paper")
		logging.info("CITATION: " +str(citation))
		handle = Entrez.efetch(db="pmc", id=citation, rettype='full', retmode="xml")
		xml_record = handle.read() #xml str
		#print(xml_record)
		logging.info("* got xml record")
		main_text, abstract_check, whole_article_check = parsePMC(xml_record, pmid)
		for yn in abstract_check:
			all_abstract_check.append(yn)
		for yn in whole_article_check:
			all_article_check.append(yn)
		logging.info("* ready to print it")
		save_path = '/home/hclent/data/'+(str(pmid))+'/' #must save to data, in proper file
		completeName = os.path.join(save_path, ('self_'+str(pmid)+'.txt')) #save self_paper as self_1234.txt
		sys.stdout = open(completeName, "w")
		print(main_text)
		i += 1
		time.sleep(3)
	logging.info("got documents: done in %0.3fs." % (time.time() - t0))
	return all_abstract_check, all_article_check


#retrieve the txt for input papers as well
def getSelfText(pmid):
	pmcid_list = [getAlternativeId(pmid)] #pmid --> list of 1 pmcid
	fake_list = ["None"]
	if pmcid_list == fake_list:
		logging.info("there wasn't a matching PMCID for the input PMID")
	else:
		logging.info(pmcid_list)
		get_self_ContentPMC(pmid, pmcid_list)  # get self text

# getSelfText("18952863")
# getSelfText("18269575")

# stuff = getMainInfo("26503504")
# print(stuff)
# pmc_ids = getCitationIDs("26503504")
# #pmc_ids = [3159747, 3122376, 3117012] #middle PMC here is bad
# pmc_titles, pmc_authors, pmc_journals, pmc_dates, pmc_urls = getCitedInfo(pmc_ids)
# all_abstract_check, all_article_check = getContentPMC("26503504", pmc_ids)
# main_info = list(zip(pmc_titles, pmc_authors, pmc_journals, pmc_dates, pmc_urls, all_abstract_check, all_article_check))
#
# logging.info(main_info)
# logging.info("############")
# for tup in main_info:
# 	title = tup[0]
# 	author = tup[1]
# 	journal = tup[2]
# 	pubdate = tup[3]
# 	url = tup[4]
# 	abstract = tup[5]
# 	article = tup[6]
# 	logging.info(title)
# 	logging.info("############")
# 	logging.info(author)
# 	logging.info("############")
# 	logging.info(journal)
# 	logging.info("############")
# 	logging.info(pubdate)
# 	logging.info("############")
# 	logging.info(url)
# 	logging.info("############")
# 	logging.info(abstract)
# 	logging.info("############")
# 	logging.info(article)
# 	logging.info("##############################")



################### Notes ##############
#Rarely, the XML will return this:
	# <?xml version="1.0"?>
	# <!DOCTYPE pmc-articleset PUBLIC "-//NLM//DTD ARTICLE SET 2.0//EN" "http://dtd.nlm.nih.gov/ncbi/pmc/articleset/nlm-articleset-2.0.dtd">
	# <pmc-articleset><Reply Id="24948109" error="The following PMCID is not available: 24948109"/></pmc-articleset>
#Thus pasrePMC() has exception handeling for this

# Beginning August, NCBI is having connection errors: "pmc is not db". Must add exception handling for it