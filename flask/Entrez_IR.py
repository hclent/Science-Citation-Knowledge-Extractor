from __future__ import print_function
from Bio import Entrez
import time, sys, pickle
from time import sleep
import datetime
import os.path, logging, json, re
import xml.etree.ElementTree as ET
from collections import defaultdict



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
	t0 = time.time()
	handle = Entrez.esummary(db="pubmed", id=pmid)
	record = Entrez.read(handle)
	title = [record[0]["Title"]] #make a list
	authors = [record[0]["AuthorList"]]
	journal = [record[0]["FullJournalName"]]
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
		handle = Entrez.esummary(db="pmc", id=citation)
		record = Entrez.read(handle)
		t = record[0]["Title"]
		a = record[0]["AuthorList"]
		j = record[0]["FullJournalName"]
		d = record[0]["PubDate"]
		u = "http://www.ncbi.nlm.nih.gov/pmc/articles/PMC"+citation
		pmc_titles.append(t)
		pmc_authors.append(a)
		pmc_journals.append(j)
		pmc_dates.append(d)
		pmc_urls.append(u)
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
		logging.DEBUG("The following PMCID is not available")
	try:
		#Get main text and add to doc
		text = root.findall('.//p')
		for t in text:
			full_text = ("".join(t.itertext()))
			main_text.append(full_text)
		#print("* Got main text")
	except Exception as e:
		logging.DEBUG("Only gave us the absract")
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


# self_info = getMainInfo(my_pmid)
# print()
#"18952863, 18269575"
#p1 = "18952863"
#pmc_ids1 = getCitationIDs(p1)
# print(pmc_ids)
# print()
# amount = len(pmc_ids)
# print("THERE ARE " + str(amount) + " DOCUMENTS")
# print()
#pmc_titles1, pmc_authors1, pmc_journals1, pmc_dates1, pmc_urls1 = list(getCitedInfo(pmc_ids1))

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

def journals_vis(journals, dates):
	years_list = []
	#print(dates) #step1 : get years
	for d in dates:
		y = re.sub('.[A-Z]{1}[a-z]{2}(.?\d{1,2})?', '', d) #delete month and day
		years_list.append(y)
	#print(years_list)

	years_range = (years_list[-1], years_list[0]) #define years range with (oldest, newest)
	#print(years_range)

	#Associate journals with years
	journal_year = list(zip(journals, years_list)) #('Scientific Reports', '2016')
	#print(journal_year)

	#Dictionary with "Journal": [year, year]
	#For looking up the years
	jyDict = defaultdict(list)
	i = 0
	for j in journals:
		if j == (journal_year[i][0]):
			jyDict[j] += [journal_year[i][1]]
			i+=1
	#print(jyDict)

	#Dictionary with "Journal": Number-of-publications
	#For looking up the total
	journalsTotalDict = defaultdict(lambda: 0)
	sum = 0
	for j in journals:
		journalsTotalDict[j] += 1
		sum +=1
	#print(journalsTotalDict)
	unique_journals = list(journalsTotalDict.keys())
	#print(unique_journals)

	publication_data = []
	for j in unique_journals:
		#print(j)
		#Initiate the dictionary for this journal
		journal_data = {
			"name": j,
			"articles": [], #[[year, number], [year, number]]
			"total": journalsTotalDict[j]   #total can get from journalsTotalDict with key (total is value)
		}
		#print("Years a paper was in this journal: "+ str(jyDict[j]))
		for year in range(int(years_range[0]), int(years_range[1]) + 1):
			#print("checking " +str(year) +" ...")
			sum = 0
			for entry in jyDict[j]:
				#print(" ... against "+str(entry))
				if year == int(entry):
					#print("The years match so I'm going to count now")
					sum+=1
				year_sum = [year, sum]
				#print(year_sum)
			journal_data["articles"].append(year_sum)

		publication_data.append(journal_data)
	publication_data = re.sub('\'', '\"', str(publication_data)) #json needs double quotes, not single quotes
	return publication_data


