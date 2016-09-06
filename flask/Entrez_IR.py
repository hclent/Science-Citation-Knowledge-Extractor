from __future__ import print_function
from Bio import Entrez
import time, sys, pickle
from time import sleep
import datetime
import os.path, logging, json, re
import xml.etree.ElementTree as ET



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
#Output: basic info on pmid as lists
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
	logging.info("self info: done in %0.3fs." % (time.time() - t0))
	return list(zip(title, authors, journal))


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
	#Get abstract and add to doc
	try:
		abstract = root.find('.//abstract')
		full_abs = ("".join(abstract.itertext()))
		#print("* Got abstract")
		main_text.append(full_abs)
	except Exception as e:
		logging.info("The following PMCID is not available")
	try:
		#Get main text and add to doc
		text = root.findall('.//p')
		for t in text:
			full_text = ("".join(t.itertext()))
			main_text.append(full_text)
		#print("* Got main text")
	except Exception as e:
		logging.info("Only gave us the absract")
	return main_text




#Input: pmid and the list of pmcids citing it
#For each citing pmc_id, this function gest the xml, which is then parsed by parsePMC()
#Output: Journal texts for each pmcid
def getContentPMC(pmid, pmcids_list):
	t0 = time.time()
	i = 1
	dirname = '/home/hclent/data/'+pmid
	try:
		os.makedirs(dirname) #creates folder named after pmid
	except OSError:
		if os.path.isdir(dirname):
			pass
		else:
			raise
	for citation in pmcids_list:
		handle = Entrez.efetch(db="pmc", id=citation, rettype='full', retmode="xml")
		xml_record = handle.read() #xml str
		#print(xml_record)
		#print("* got xml record")
		main_text = parsePMC(xml_record, pmid)
		#print("* ready to print it")
		save_path = '/home/hclent/data/'+(str(pmid))+'/' #must save to data, in proper file
		completeName = os.path.join(save_path, (str(pmid)+'_'+str(i)+'.txt'))
		sys.stdout = open(completeName, "w")
		print(main_text)
		i += 1
		time.sleep(3)
	logging.info("got documents: done in %0.3fs." % (time.time() - t0))


# stuff = getMainInfo("22555442")
# print(stuff)
# pmc_ids = getCitationIDs("22555442")
# print(pmc_ids)
# pmc_titles, pmc_authors, pmc_journals, pmc_dates, pmc_urls = getCitedInfo(pmc_ids)
# print(pmc_titles)

################### Notes ##############
#Rarely, the XML will return this:
	# <?xml version="1.0"?>
	# <!DOCTYPE pmc-articleset PUBLIC "-//NLM//DTD ARTICLE SET 2.0//EN" "http://dtd.nlm.nih.gov/ncbi/pmc/articleset/nlm-articleset-2.0.dtd">
	# <pmc-articleset><Reply Id="24948109" error="The following PMCID is not available: 24948109"/></pmc-articleset>
#Thus pasrePMC() has exception handeling for this

# Beginning August, NCBI is having connection errors: "pmc is not db". Must add exception handling for it