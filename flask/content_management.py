from processors import * #pyProcessors
import os.path, time, re, logging
from Entrez_IR import * #mine
from multi_preprocess import * #mine
from lsa1 import * #mine
from lda1 import * #mine
from journalvis import * #mine


## Supporting functions for app.py


#Create log
logging.basicConfig(filename='.app.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


#Set a PROCESSORS_SERVER environment variable.
#It may take a minute or so to load the large model files.
def connect_to_Processors(port_num):
  path = '/home/hclent/anaconda3/envs/py34/lib/python3.4/site-packages/processors/processors-server.jar'
  api = ProcessorsAPI(port=port_num, jar_path=path, keep_alive=True, jvm_mem="-Xmx16G")
  logging.info('Connected to pyProcessors')
  #Initialize the bionlp annotator by passing it a doc
  init_doc = api.bionlp.annotate("The mitochondria is the powerhouse of the cell.")
  return api


#User can enter in as many pubmed ids as they want into text box
#This method creates a list of them
def multiple_pmid_input(user_input):
	logging.info('cleaning user input')
	clean = re.sub('\,', ' ', user_input)
	ids = clean.split() #list of pmids
	return ids


#If pmid (user input) in the database, just get main_info (authors, journals, ect)
def run_IR_in_db(user_input):
	logging.info('PMID is in the database')
	self_info = getMainInfo(user_input)
	pmc_ids = getCitationIDs(user_input)
	target_title, target_authors, target_journals, target_dates, target_urls = getCitedInfo(pmc_ids)
	main_info = list(zip(target_title, target_authors, target_urls))
	return main_info, target_journals, target_dates


#If pmid (user input) NOT in the db, get main_info AND scrape XML for abstracts and texts
def run_IR_not_db(user_input):
	logging.info('PMID is NOT in the database')
	self_info = getMainInfo(user_input)
	pmc_ids = getCitationIDs(user_input)
	target_title, target_authors, target_journals, target_dates, target_urls = getCitedInfo(pmc_ids)
	main_info = list(zip(target_title, target_authors, target_urls))
	#Get XML
	getContentPMC(user_input, pmc_ids)
	return main_info, target_journals, target_dates


def print_journalvis(journals, dates, user_input):
	publication_data = journals_vis(journals, dates)
	print(publication_data)
	logging.info('Printing JOURNALS to json')
	save_path = '/home/hclent/data/'+str(user_input)+'/'
	completeName = os.path.join(save_path, ('journals_'+(str(user_input))+'.json'))
	with open(completeName, 'w') as outfile:
		json.dump(publication_data, outfile)



#Take pmid_n.txt and get an annotated document, as well as lemmas and named entities
#This method is for user_input NOT already in DB, need to make json, for in DB, no need to make JSON
def do_ALL_multi_preprocessing(user_input):
	logging.info('Beginning multiprocessing for NEW docs')
	t1 = time.time()
	docs = retrieveDocs(user_input)
	multiprocess(docs)
	biodocs = retrieveBioDocs(user_input)
	data_samples, nes_list = loadBioDoc(biodocs)
	logging.info("Execute everything: done in %0.3fs." % (time.time() - t1))
	return data_samples, nes_list


#Take annotated docs and return data and nes
#This method is for user_input that IS already in the DB
def do_SOME_multi_preprocessing(user_input):
	logging.info('Beginning multiprocessing for PRE-EXISTING docs')
	#t1 = time.time()
	biodocs = retrieveBioDocs(user_input)
	data_samples, nes_list = loadBioDoc(biodocs)
	#logging.info("Execute everything: done in %0.3fs." % (time.time() - t1))
	return data_samples, nes_list



def run_lsa1(user_input, data_samples, k):
	logging.info('Beginning Latent Semantic Analysis')
	tfidf, tfidf_vectorizer = get_tfidf(data_samples)
	jsonDict = do_LSA(tfidf, tfidf_vectorizer, k) #need to make this an option
	return jsonDict

def print_lsa(user_input, jsonDict):
	#Save the json for @app.route('/reslsa/')
	logging.info('Printing LSA to json')
	save_path = '/home/hclent/data/'+str(user_input)+'/'
	completeName = os.path.join(save_path, ('lsa_'+(str(user_input))+'.json'))
	with open(completeName, 'w') as outfile:
		json.dump(jsonDict, outfile)


def run_lda1(data_samples): #set at defulat k=3, number of words=6
	logging.info('Beginning Latent Dirichlet Allocation')
	tfidf, tfidf_vectorizer = get_tfidf(data_samples)
	lda = fit_lda(tfidf)
	jsonLDA = topics_lda(tfidf_vectorizer, lda)
	#print(jsonDict)
	return jsonLDA


def print_lda(user_input, jsonLDA):
	#Save the json for @app.route('/reslda/')
	logging.info('Saving LDA to JSON')
	save_path = '/home/hclent/data/'+str(user_input)+'/'
	completeName = os.path.join(save_path, ('lda_'+(str(user_input))+'.json'))
	with open(completeName, 'w') as outfile:
		json.dump(jsonLDA, outfile)