from processors import * #pyProcessors
from Entrez_IR import * #mine
from multi_preprocess import * #mine
from lsa1 import * #mine
import os.path, time, re


#Set a PROCESSORS_SERVER environment variable.
#It may take a minute or so to load the large model files.
def connect_to_Processors(port_num):
  path = '/home/hclent/anaconda3/envs/pyProcessors/lib/python3.4/site-packages/py34/processors-server.jar'
  api = ProcessorsAPI(port=port_num, jar_path=path, keep_alive=True, jvm_mem="-Xmx16G")
  logging.info('Connected to pyProcessors')
  #Initialize the bionlp annotator by passing it a doc
  init_doc = api.bionlp.annotate("The mitochondria is the powerhouse of the cell.")
  return api


#User can enter in as many pubmed ids as they want into text box
#This method creates a list of them
def multiple_pmid_input(user_input):
    clean = re.sub('\,', ' ', user_input)
    ids = clean.split() #list of pmids
    return ids


#If pmid (user input) in the database, just get main_info (authors, journals, ect)
def run_IR_in_db(user_input):
	self_info = getMainInfo(user_input)
	pmc_ids = getCitationIDs(user_input)
	target_title, target_authors, target_journals, target_urls = getCitedInfo(pmc_ids)
	main_info = list(zip(target_title, target_authors, target_urls))
	return main_info, target_journals


#If pmid (user input) NOT in the db, get main_info AND scrape XML for abstracts and texts
def run_IR_not_db(user_input):
	self_info = getMainInfo(user_input)
	pmc_ids = getCitationIDs(user_input)
	target_title, target_authors, target_journals, target_urls = getCitedInfo(pmc_ids)
	main_info = list(zip(target_title, target_authors, target_urls))
	#Get XML
	getContentPMC(user_input, pmc_ids)
	return main_info, target_journals


# #Take pmid_n.txts and get stuff!
# #Update: for not in DB, need to make json, for in DB, no need to make JSON
# def do_preprocessing(num, user_input, api):
# 	loadDocuments(num, user_input, api) #loads txt, pre-processes, and dumps to JSON
# 	print(" * Converted txt to Docs ...")
# 	data_samples, ner_list = loadBioDoc(num, user_input) #loads json
# 	print(" * Got info from Docs ...")
# 	return data_samples, ner_list
#
#
# def already_have_preproc(num, user_input):
# 	data_samples, ner_list = loadBioDoc(num, user_input) #loads json
# 	print("* Got info from Docs ...")
# 	return data_samples, ner_list


#Annotate docs to json, and return data and nes
#This method is for user_input NOT already in DB, need to make json, for in DB, no need to make JSON
def do_ALL_multi_preprocessing(user_input, api):
	t1 = time.time()
	docs = retrieveDocs(user_input)
	multiprocess(docs)
	biodocs = retrieveBioDocs(user_input)
	data_samples, nes_list = loadBioDoc(biodocs)
	print("Execute everything: done in %0.3fs." % (time.time() - t1))
	return data_samples, nes_list


#Take annotated docs and return data and nes
#This method is for user_input that IS already in the DB
def do_SOME_multi_preprocessing(user_input, api):
	t1 = time.time()
	biodocs = retrieveBioDocs(user_input)
	data_samples, nes_list = loadBioDoc(biodocs)
	print("Execute everything: done in %0.3fs." % (time.time() - t1))
	return data_samples, nes_list



def run_lsa1(user_input, data_samples, k):
	tfidf, tfidf_vectorizer = get_tfidf(data_samples)
	jsonDict = do_LSA(tfidf, tfidf_vectorizer, k) #need to make this an option
	print(jsonDict)
	print(" * Generated json for LSA visualization !")
	# Unsure whether or not I should save JSON or just pass it straight to Flask ... 
	# save_path = '/Users/hclent/Desktop/webdev-biotool/flask/static/' #must save to static
	# completeName = os.path.join(save_path, ('vis_'+(str(user_input))+'.json'))
	# with open(completeName, 'w') as outfile:
	# 	json.dump(jsonDict, outfile, default=dumper, indent=2)
	return jsonDict