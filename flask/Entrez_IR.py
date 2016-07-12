from __future__ import print_function
from Bio import Entrez
import time, sys, pickle
from time import sleep
import datetime
import os.path
import xml.etree.ElementTree as ET


#Entrez Information Retrieval
#This code uses BioPython to access NCBI's API (Entrez)
#From the NCBI API, this code references PubMed and PubMedCentral
#With an input PubMed ID's (pmids), the code will retrieve information about the original pmid,
#and information about pubplications that cite this pmid via PubMedCentral ID's (pmcids)


Entrez.email = "hclent1@gmail.com" 
Entrez.tool = "MyInfoRetrieval"

my_pmid = 18952863

#Input: pmid
#Output: basic info on pmid as lists
def getMainInfo(pmid):
	t0 = time.time()
	handle = Entrez.esummary(db="pubmed", id=pmid)
	record = Entrez.read(handle)
	title = [record[0]["Title"]] #make a list
	authors = [record[0]["AuthorList"]]
	journal = [record[0]["FullJournalName"]]
	print("self info: done in %0.3fs." % (time.time() - t0))
	return list(zip(title, authors, journal))




#Input: Pmid
#Output: list of pmcids of the articles citing you
def getCitationIDs(pmid): #about the same speed as MainCrawl.py
	t0 = time.time()
	results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", from_uid=pmid))
	pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
	return pmc_ids #list
	print("get Citation PMCIDs: done in %0.3fs." % (time.time() - t0))
	time.sleep(3)



#Input: Citing pmcids
#Output: Basic info about these pmcids
def getCitedInfo(pmcid_list): 
	t0 = time.time()
	pmc_titles = []
	pmc_authors = []
	pmc_journals = []
	pmc_urls = []
	i = 1
	for citation in pmcid_list:
		print("citation no. " + str(i) + " ...")
		handle = Entrez.esummary(db="pmc", id=citation)
		record = Entrez.read(handle)
		t = record[0]["Title"]
		a = record[0]["AuthorList"]
		j = record[0]["FullJournalName"]
		u = "http://www.ncbi.nlm.nih.gov/pmc/articles/PMC"+citation
		pmc_titles.append(t)
		pmc_authors.append(a)
		pmc_journals.append(j)
		pmc_urls.append(u)
		time.sleep(3)
		i += 1
	#main_info = (list(zip(pmc_titles, pmc_authors, pmc_journals, pmc_urls)))
	print("get citations info: done in %0.3fs." % (time.time() - t0))
	return pmc_titles, pmc_authors, pmc_journals, pmc_urls


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
		print("The following PMCID is not available")
	try:
		#Get main text and add to doc
		text = root.findall('.//p')
		for t in text:
			full_text = ("".join(t.itertext()))
			main_text.append(full_text)
		#print("* Got main text")
	except Exception as e:
		print("Only gave us the absract")
	return main_text




#Input: pmid and the list of pmcids citing it
#For each citing pmc_id, this function gest the xml, which is then parsed by parsePMC()
#Output: Journal texts for each pmcid
def getContentPMC(pmid, pmcids_list):
	t0 = time.time()
	i = 1
	for citation in pmcids_list:
		handle = Entrez.efetch(db="pmc", id=citation, rettype='full', retmode="xml")
		xml_record = handle.read() #xml str
		#print(xml_record)
		#print("* got xml record")
		main_text = parsePMC(xml_record, pmid)
		#print("* ready to print it")
		save_path = '/Users/hclent/Desktop/webdev-biotool/flask/data/' #must save to data
		completeName = os.path.join(save_path, (str(pmid)+'_'+str(i)+'.txt'))
		sys.stdout = open(completeName, "w")
		print(main_text)
		i += 1
		time.sleep(3)
	print("got documents: done in %0.3fs." % (time.time() - t0))




# self_info = getMainInfo(my_pmid)
# print()
# pmc_ids = getCitationIDs(my_pmid)
# print(pmc_ids)
# print()
# amount = len(pmc_ids)
# print("THERE ARE " + str(amount) + " DOCUMENTS")
# print()
# main_info = list(getCitedInfo(pmc_ids))
# print(main_info)
# with open('p18952863.pickle', 'wb') as f:
# 	pickle.dump(main_info, f)
# print("its been pickled yo!")
# journals = []
# with open('p18952863.pickle', 'rb')as f:
# 	data = pickle.load(f)
# for citations in data:
# 	j = citations[2]
# 	journals.append(j)
# print(journals)
# 	print("PICKLED VER YO")
# 	print(data)
#getContentPMC(my_pmid, pmc_ids)

################### Notes ##############
#Rarely, the XML will return this:
	# <?xml version="1.0"?>
	# <!DOCTYPE pmc-articleset PUBLIC "-//NLM//DTD ARTICLE SET 2.0//EN" "http://dtd.nlm.nih.gov/ncbi/pmc/articleset/nlm-articleset-2.0.dtd">
	# <pmc-articleset><Reply Id="24948109" error="The following PMCID is not available: 24948109"/></pmc-articleset>
#Thus pasrePMC() has exception handeling for this

