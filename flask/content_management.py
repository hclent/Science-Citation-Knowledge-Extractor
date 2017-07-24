from processors import * #pyProcessors
from flask import Flask
import os.path, time, re, logging, pickle, json, codecs, arrow
import operator
from configapp import app, engine, connection, inputPapers, citations, queries
from database_management import * #mine
from Entrez_IR import * #mine
from multi_preprocess import * #mine
from lsa1 import * #mine
from lda1 import * #mine
from journalvis import * #mine
from nes import * #mine
from kmeans1 import * #mine
from naive_cosineSim import * #mine
from fasttext import * #mine
from fgraph import * #mine
from fgraph2json import embedding_json #mine
from cache_lemma_nes import load_lemma_cache #mine


############## CONFIG ######################################

################## LOGGING #######################################################

logging.basicConfig(filename=app.config['PATH_TO_LOG'],level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logging.debug("Why dost thine tail -F break?")

############# PROCESSORS SERVER ##################################################

#It may take a minute or so to load the large model files.
def connect_to_Processors(port_num):
  path = (app.config['PATH_TO_JAR'])
  api = ProcessorsAPI(port=port_num, jar_path=path, keep_alive=True, jvm_mem="-Xmx100G")
  logging.info('Connected to pyProcessors')
  #Initialize the bionlp annotator by passing it a doc
  init_doc = api.bionlp.annotate("The mitochondria is the powerhouse of the cell.")
  return api

################ WORD EMBEDDING MODEL ############## (Global var)

#This is really slow to load :/ Maybe shouldn't be a global var?
# Don't want to re-load for each analysis though...
fasttext_model = load_model(app.config['PATH_TO_FASTTEXT_MODEL'])

################### INPUT #########################################################

#User can enter in as many pubmed ids as they want into text box
#This method creates a list of them
def multiple_pmid_input(user_input):
	logging.info('cleaning user input')
	clean = re.sub('\,', ' ', user_input)
	ids = clean.split() #list of pmids
	return ids

def flatten(listOfLists):
    return list(chain.from_iterable(listOfLists))

################### DATABASE #####################################################
#If the pmid is NOT in the db, we also need to getSelfText and write that info to db
def scrape_and_write_Input(user_input, conn):
	logging.info('retrieve Self text, and write to db')
	self_pmcid, self_abstract_check, self_article_check = getSelfText(user_input) #Entrez_IR function
	logging.info("PMICD: ")
	logging.info(self_pmcid)
	updateInputPapers(user_input, self_pmcid, self_abstract_check, self_article_check, conn)  # put getSelfText into database

#input: user_input pmid
#output: list of dicts with annotation checks [{"pmcid": 123, "annotated": yes}]
def annotation_check(user_input, conn):
	a_check = []
	pmc_ids = db_citation_pmc_ids(user_input, conn) #Used to use getCitationIDs(user_input) here but updated to using my own db
	for citation in pmc_ids:
		#print(citation)
		annotationDict = {"pmcid": citation, "annotated": []}

		prefix = str(citation[0:3])  # folder for first 3 digits of pmcid
		suffix = prefix + '/' + str(citation[3:6])  # folder for second 3 digits of pmcid nested in prefix
		filename = suffix + '/' + str(citation) + '.json'
		full_filename = os.path.join((app.config['PATH_TO_CACHE']) ,filename)
		with open(full_filename) as data_file:
			data = json.load(data_file)
			if "error annotating document" in data["text"][:25]:
				annotationDict["annotated"].append("no")
			else:
				annotationDict["annotated"].append("yes")

		a_check.append(annotationDict)
	return a_check

#Used by both IR_not_in_db and IR_in_db to 1) add new pmcids to cited db and 2) duplicate entries when necessary
#Input is the "citation" (dict) result of the result "allCitationsInfo = getCitedInfo(pmc_ids, user_input)" in for loop
#Updated to sqlalchemy
def new_or_copy_db(citation, conn): #citation is a dict
	if "annotated" in citation:
		logging.info("this entry has already been annotated before. just copy.")
		date = str(arrow.now().format('YYYY-MM-DD'))
		pmcid = citation["pmcid"]
		title = str(citation["pmc_titles"][0])
		author = str(citation["pmc_authors"][0])
		journal = str(citation["pmc_journals"][0])
		pubdate = str(citation["pmc_dates"][0])
		citesPmid = str(citation["citesPmid"])
		url = str(citation["pmc_urls"][0])
		abstract = str(citation["abstract_check"][0])
		whole_article = str(citation["article_check"][0])
		sents = str(citation["sents"][0])
		tokens = str(citation["tokens"][0])
		annotated = str(citation["annotated"][0])

		update = citations.insert().\
			values(dict(datestamp = date, pmcid=pmcid, title=title, author=author, journal=journal, pubdate=pubdate,
						citesPmid=citesPmid, url=url, abstract=abstract, whole_article=whole_article, sents=sents,
						tokens=tokens, annotated=annotated ))
		conn.execute(update)

	if "annotated" not in citation:
		logging.info("this entry is brand new, never annotated")
		date = str(arrow.now().format('YYYY-MM-DD'))
		pmcid = str(citation["pmcid"])
		title = str(citation["pmc_titles"][0])
		authorlist = citation["pmc_authors"][0]
		s = ', '
		author = str(s.join(authorlist))
		journal = str(citation["pmc_journals"][0])
		pubdate = str(citation["pmc_dates"][0])
		citesPmid = str(citation["citesPmid"])
		url = str(citation["pmc_urls"][0])

		update = citations.insert().\
			values(dict(datestamp = date, pmcid=pmcid, title=title, author=author, journal=journal, pubdate=pubdate,
						citesPmid=citesPmid, url=url))
		conn.execute(update)


#If pmid (user input) NOT in the db, get main_info AND scrape XML for abstracts and texts
#Write self_info to inputPmids db
#Write allCitationsInfo to citations db
#Update citations db with abstract_check and whole_article_check
#Doesn't re-retrieve information for citations previously scraped (does update db if necessary)
# Updated to sqlalchemy
def run_IR_not_db(user_input, conn):
	logging.info('PMID is NOT in the inputPapers database')
	self_info = getMainInfo(user_input)
	pmc_ids = getCitationIDs(user_input)
	num_citations = len(pmc_ids)
	logging.info("Writing self_info to inputPapers db")
	#write self_info to "inputPapers" db
	for tup in self_info:
		title = tup[0]
		s = ', '
		author = str(s.join(tup[1]))
		journal = tup[2]
		pubdate = tup[3]
		url = tup[4]
		date = str(arrow.now().format('YYYY-MM-DD'))

		update = inputPapers.insert().\
			values(dict(datestamp=date, pmid=user_input, title=title, author=author, journal=journal, pubdate=pubdate,
						url=url, num_citations=num_citations))
		conn.execute(update)

	#Retrieve the input paper if avaliable and update db
	scrape_and_write_Input(user_input, conn)

	#Now retrieve citations
	logging.info("Get basic info about the citations")
	# Previously unseen pmcids only in allCitationsInfo.
	# Previously seen pmcids are copied to db for new pmid in getCitedInfo
	allCitationsInfo = getCitedInfo(pmc_ids, user_input, conn) #output: list of dictionaries [{pmid: 1234, author: human, ...}]
	logging.info("Write basic citation info to citations db")
	for citation in allCitationsInfo:
		logging.info(citation)
		new_or_copy_db(citation, conn)

	#Get content and update citations db
	contentDictList = getContentPMC(pmc_ids, user_input, conn)
	for citation in contentDictList:
		pmcid = str(citation["pmcid"])
		citesPmid = str(citation["citesPmid"])
		abstract = str(citation["all_abstract_check"][0])
		whole_article = str(citation["all_article_check"][0])

		up = citations.update().\
			where(citations.c.pmcid == pmcid).\
			where(citations.c.citesPmid == citesPmid).\
			values(dict(abstract=abstract, whole_article=whole_article))
		conn.execute(up)


# If pmid (user input) in the inputPapers database, check for new papers
# If ther are new citing papers, get those and update the db appropriately
# Else if there are no new papers, don't need to do anything.
# Updated to sqlalchemy
def run_IR_in_db(user_input, conn):
	logging.info('PMID is in the database')
	# Check for new papers:
	num_in_db = db_input_citations_count(user_input, conn) #checks MY db
	pmc_ids = getCitationIDs(user_input) #checks ENTREZ DB
	num_current = len(pmc_ids)
	#If there are new papers,
	if int(num_current) > int(num_in_db): #TODO change this back to > after i've fixed authors problem
		need_to_annotate = 'yes'
		#print("there are new citations!", (num_current, num_in_db))
		logging.info("num_in_db: ", num_in_db)
		logging.info("num_current: ", num_current)
		#update number of citations in inputPaper db
		update = inputPapers.update().\
			where(inputPapers.c.pmid == user_input).\
			values(num_citations=num_current)
		conn.execute(update)

		#now get the new citation info
		allCitationsInfo = getCitedInfo(pmc_ids, user_input, conn)  # output: list of dictionaries [{pmid: 1234, author: human, ...}] #skips duplicates
		logging.info("Write basic citation info to citations db for new papers")
		for citation in allCitationsInfo:
			new_or_copy_db(citation, conn)

		#Get content and update citations db
		logging.info("now get the content for the new stuff")
		contentDictList = getContentPMC(pmc_ids, user_input, conn)
		for citation in contentDictList:
			pmcid = str(citation["pmcid"])
			citesPmid = str(citation["citesPmid"])
			abstract = str(citation["all_abstract_check"][0])
			whole_article = str(citation["all_article_check"][0])

			up = citations.update().\
				where(citations.c.pmcid == pmcid).\
				where(citations.c.citesPmid == citesPmid).\
				values(dict(abstract=abstract, whole_article=whole_article))
			conn.execute(up)

	else:
		logging.info("no new papers, nothing to do here folks")
		need_to_annotate = 'no'
		pass
	return need_to_annotate


#Data for populating statistics page in app
def get_statistics(query, conn):
	total_pubs, unique_pubs, abstracts, whole, sentences, words =  db_query_statistics(query, conn)
	statistics = [total_pubs, unique_pubs, abstracts, whole, sentences, words]
	return statistics


#use query to get info about input papers
def statsSelfInfo(query, conn):
    input_click_citations = []
    pmid_list = query.split('+')  # list of string pmids
    for user_input in pmid_list:
        apa = db_inputPapers_retrieval(user_input, conn)
        url = "https://www.ncbi.nlm.nih.gov/pubmed/"+str(user_input)
        href_label = (apa, url) #store apa and url as a tuple
        input_click_citations.append(href_label) #then append to list
    return input_click_citations


#take a query and generate x and y datapoints for pubs x year bar chart in "Stats" tab
#I'm sorry this is so hacky and ugly :'(
def stats_barchart(query, conn):
	logging.info("STATISTICS: stacked bar chart initializing ... ")
	pmid_list = query.split('+') #list of string pmids

	#maximum of 5 stacked barcharts
	x0 = []
	x1 = []
	x2 = []
	x3 = []
	x4 = []

	y0 = []
	y1 = []
	y2 = []
	y3 = []
	y4 = []

	n0 = []
	n1 = []
	n2 = []
	n3 = []
	n4 = []

	for user_input in pmid_list:
		journals = []
		dates = []
		db_journals, db_dates = db_bar_chart(user_input, conn)
		for j in db_journals:
			journals.append(j)
		for d in db_dates:
			dates.append(d)
		x, y = statistics_dates_barchart(journals, dates, query, conn)
		logging.info(x)
		logging.info(y)

		#should at least have ONE paper...
		#I'm SO sorry world for this terribly terribly hacky sollution :'(
		if user_input == pmid_list[0]:
			x0.append(x)
			y0.append(y)
			n0.append(user_input)
		try:
			if user_input == pmid_list[1]:
				x1.append(x)
				y1.append(y)
				n1.append(user_input)
		except Exception as p1:
			pass
		try:
			if user_input == pmid_list[2]:
				x2.append(x)
				y2.append(y)
				n2.append(user_input)
		except Exception as p2:
			pass
		try:
			if user_input == pmid_list[3]:
				x3.append(x)
				y3.append(y)
				n3.append(user_input)
		except Exception as p3:
			pass
		try:
			if user_input == pmid_list[4]:
				x4.append(x)
				y4.append(y)
				n4.append(user_input)
		except Exception as p4:
			pass

	x0 = flatten(x0)
	x1 = flatten(x1)
	x2 = flatten(x2)
	x3 = flatten(x3)
	x4 = flatten(x4)
	y0 = flatten(y0)
	y1 = flatten(y1)
	y2 = flatten(y2)
	y3 = flatten(y3)
	y4 = flatten(y4)
	n0 = str(''.join(n0))
	n1 = str(''.join(n1))
	n2 = str(''.join(n2))
	n3 = str(''.join(n3))
	n4 = str(''.join(n4))

	return x0, x1, x2, x3, x4, y0, y1, y2, y3, y4, n0, n1, n2, n3, n4


############ PROCESSING DOCS --> BIODOCS ############################################
#Take pmcid.txt and get an annotated document, as well as lemmas and named entities
#Doesn't re-annotated documents that have already been annotated.
#Updated to sqlalchemy
#This is only called when there are new documents to annotate :)
def do_multi_preprocessing(user_input, conn):
	logging.info('Beginning multiprocessing for NEW (unprocessed) docs')
	t1 = time.time()
	docs = retrieveDocs(user_input, conn)
	multiprocess(docs) #if docs is empty [], this function just passes :)

	#Now update annotated_check
	a_check = annotation_check(user_input, conn)

	for a in a_check: #{"pmcid": pmcid, "annotated": ['yes']}
		logging.info("updating the annotation checks in the db")
		pmcid = str(a["pmcid"])
		annotated = str(a["annotated"][0])

		update = citations.update().\
			where(citations.c.pmcid == pmcid).\
			where(citations.c.citesPmid == user_input).\
			values(annotated=annotated)
		conn.execute(update)

	#Now extract information from annotated documents
	biodocs = retrieveBioDocs(user_input, conn)
	biodoc_data = loadBioDoc(biodocs) #list of dictionaries[{pmid, lemmas, nes, sent_count, token_count}]
	#No problem getting biodocs or biodoc_data ... problem comes with updating db...
	#update db with sents and tokens
	for b in biodoc_data:
		update_annotations(b, user_input, conn)
	logging.info("Execute everything: done in %0.3fs." % (time.time() - t1))
	return biodoc_data

############ DATA VISUALIZATIONS #################################################

#This code is called in the function below (print_journalvis)
#Basically it forces an update of the journals vis & update to "queries" table of db.
#We would want to force an update of the journals vis if there are new papers to a previously seen query
def force_update_journals(query, conn):
	years_range = get_years_range(query, conn)  # need range for ALL journals, not just last one
	logging.info(years_range)
	publication_data, range_info = journals_vis(years_range, query, conn)  # range info = [('2008', '2016'), 165, 48]
	logging.info(range_info)
	journal_years = range_info[0]
	logging.info(journal_years)
	q = '+'
	logging.info(q)
	range_years = str(q.join(journal_years))
	logging.info("range years: " + range_years)
	unique_publications = range_info[1]
	unique_journals = range_info[2]
	logging.info(range_info)

	logging.info('Printing JOURNALS to JSON')
	pmid_list = query.split('+')  # list of string pmids
	pmid = pmid_list[0]  # get the first
	prefix = pmid[0:3]
	suffix = pmid[3:6]

	try:
		os.makedirs(os.path.join((app.config['PATH_TO_JOURNALS']), prefix))
	except OSError:
		if os.path.isdir(os.path.join((app.config['PATH_TO_JOURNALS']), prefix)):
			pass
		else:
			raise

	try:
		os.makedirs(os.path.join((app.config['PATH_TO_JOURNALS']), prefix, suffix))
	except OSError:
		if os.path.isdir(os.path.join((app.config['PATH_TO_JOURNALS']), prefix, suffix)):
			pass
		else:
			raise

	filename = str(prefix) + '/' + str(suffix) + '/' + "journals_" + str(query) + ".json"
	save_path = (app.config['PATH_TO_JOURNALS'])

	completeName = os.path.join(save_path, filename)
	with open(completeName, 'w') as outfile:
		json.dump(publication_data, outfile)
	date = str(arrow.now().format('YYYY-MM-DD'))

	update = queries.insert().\
		values(dict(datestamp=date, query=query, range_years=range_years, unique_pubs=unique_publications,
					unique_journals=unique_journals))
	conn.execute(update)
	return range_years, unique_publications, unique_journals

#Function called for actually making the journals vis.
def print_journalvis(query, needed_to_annotate_check, conn):
	#update the cache!
	if 'yes' in needed_to_annotate_check: #update the DB ()
		logging.info("Journals: NEW docs, so need to force update")
		range_years, unique_publications, unique_journals = force_update_journals(query, conn)
	# #if nothing was annotated, check the existing file. If there isn't one, pass
	if 'yes' not in needed_to_annotate_check: ## check the record.
		#check for file and get stuff from db
		record = checkForQuery(query, conn)
		if record == 'yes':  # if its in the db, just get the important things from the db!!
			logging.info("Journals: NO new docs, retrive from db")
			range_years, unique_publications, unique_journals = getJournalsVis(query, conn)
		if record == 'empty':
			logging.info("Journals: QUERY not in db table queries! Force update.")
			range_years, unique_publications, unique_journals = force_update_journals(query, conn)
	return range_years, unique_publications, unique_journals


def vis_wordcloud(neslist, nes_categories, w_number):
	nesDict = frequency_dict(neslist, nes_categories)
	wcl = wordcloud(nesDict, int(w_number))
	return wcl

def vis_heatmapTitles(lemma_samples, years, conn):
	titles = []  # want citations instead of titles

	zipped = zip(years, lemma_samples)
	sorted_data = list(sorted(zipped, key=lambda z: int(z[0])))
	pmcids = [l[1][0] for l in sorted_data]

	for id in pmcids:
		hyperlink = db_citations_hyperlink_retrieval(id, conn)
		titles.append(hyperlink[0])
	return titles

def vis_heatmap(lemma_samples, nes_samples, nes_categories, w_number, conn):
	neslist = [n[1] for n in nes_samples]
	nesDict = frequency_dict(neslist, nes_categories)
	#everything is sorted by year inside doHeatmap
	x_docs, y_words, z_counts, years = doHeatmap(nesDict, w_number, lemma_samples, conn)
	#sorted by year in vis_HeatmapTitles the sme way as doHeatmap
	titles = vis_heatmapTitles(lemma_samples, years, conn)
	return x_docs, y_words, z_counts, titles

#TODO: vis_heatmap and vis_clustermap bascially work the same.
#TODO: share the x_docs, y_words, z_counts between the two, for faster work :)
#TODO: fix papers axis so its not smooshed together and ugly
#TODO: fix the way papers are saved here, so we have the nested directories for faster lookup
#TODO: check for clustermap file
#TODO: update clustermap cache if new citations to a query
def vis_clustermap(lemma_samples, nes_samples, nes_categories, w_number, query, conn):
	logging.info("starting clustermap")
	x, y, z, titles = vis_heatmap(lemma_samples, nes_samples, nes_categories, w_number, conn)
	logging.info("making clustermap data")
	seaData = make_seaborn_data(x, y, z)
	logging.info("saving clustermap png")
	saveName = makeClusterMap(seaData, query)
	return saveName #return filename

def vis_kmeans(lemma_samples, num_clusters, conn):
	#use query to get titles
	titles = [] #want citations instead of titles

	pmcids = [l[0] for l in lemma_samples]
	for id in pmcids:
		hyperlink = db_citations_hyperlink_retrieval(id, conn)
		titles.append(hyperlink[0])

	# ignore the pmcid's in l[0], ignore tags in l[2] and just grab the lemmas
	lemmas_for_kmeans = [l[1] for l in lemma_samples]
	hX, hasher = get_hashing(lemmas_for_kmeans)
	clusters = do_kemeans(hX, int(num_clusters)) #list of cluster assignments
	coordinates = do_NMF(hX) #dimensionality reduction for visualization
	#zip coordinates and titles
	zipped_coordinates = zip(coordinates, titles)
	x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates, titles0, titles1, titles2, titles3, titles4 = plotKmeans(zipped_coordinates, clusters) #format for Plotly scatterplot
	return x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates, titles0, titles1, titles2, titles3, titles4

#scifi visualization
#are query papers eligible to be loaded as a corpus?
#step 1: does the pmid have a pmcid? If so, it was probably scraped
#step 2: check to make sure that file actually exists
#Get the data for it
#NB: This code assumes if pmid has pmcid, then it has been scraped
def inputEligible(query, conn):
	papers = []
	values = ['paper1', 'paper2', 'paper3', 'paper4', 'paper5']
	path_to_paper = []
	pmid_list = query.split('+')  # list of string pmids
	display_title = [] #title for the dropdown menu

	for pmid in pmid_list:
	#We are assuming that if there is a pmcid for the pmid, that means it was in PubmedCentral and we scraped it
		pmcid = pmid2pmcid(pmid, conn)
		if pmcid != "NA":
			#print(pmcid)
			#print(pmid + " = " + pmcid)
			#get the pmcid of the pmid

			prefix = pmcid[0:3]
			suffix = pmcid[3:6]
			path_to_txt = str(prefix) + '/' + str(suffix) + '/' + str(pmcid) + '.txt' # look in folder that matches pmcid
			filename = os.path.join((app.config['PATH_TO_CACHE']) , path_to_txt)
			truth_value = os.path.isfile(filename)
			if truth_value is True:
				#print(filename) #NA won't exist
				papers.append(pmid)
				path_to_paper.append(filename)
				display = db_pmid_axis_label(pmid, conn)
				keep_display = display[0]
				display_title.append(keep_display)

	eligible_papers = list(zip(values, papers, path_to_paper, display_title))
	return eligible_papers

#visualization for scifi div
def vis_scifi(corpus, query, eligible_papers, conn):
	corpus_vec, color = load_corpus(corpus, eligible_papers)
	eligible_cosines = get_cosine_eligible(corpus_vec, eligible_papers)
	logging.info("eligible_cosines: done")
	data_vecs_list, pmcids_list  = load_datasamples(query)
	logging.info("data_vecs_list: done")
	cosine_list = get_cosine_list(corpus_vec, data_vecs_list)
	logging.info("cosine_list: done")
	#cosine list is the cosine sim score for each document
	sorted_combos = add_urls(cosine_list, color, pmcids_list, conn)
	logging.info("sorted combos: done")
	all_sorted_combos = add_eligible_cosines(sorted_combos, eligible_papers, eligible_cosines, conn)
	logging.info("all_sorted_combos: done")
	x, y, names, color_list = prepare_for_histogram(all_sorted_combos)
	return x, y, names, color_list

############ TOPIC MODELING ############################################

#### L S A #####
def run_lsa1(lsa_lemmas, k):
	logging.info('Beginning Latent Semantic Analysis')
	tfidf, tfidf_vectorizer = get_tfidf(lsa_lemmas)
	jsonDict = do_LSA(tfidf, tfidf_vectorizer, k) #need to make this an option
	return jsonDict

def print_lsa(query, jsonDict, k):
	logging.info('Printing LSA to JSON')
	save_path = (app.config['PATH_TO_LSA'])

	#create folders if don't exist
	pmid_list = query.split('+')  # list of string pmids
	pmid = pmid_list[0]  # get the first
	prefix = pmid[0:3]
	suffix = pmid[3:6]

	try:
		os.makedirs(os.path.join((app.config['PATH_TO_LSA']), prefix))
	except OSError:
		if os.path.isdir(os.path.join((app.config['PATH_TO_LSA']), prefix)):
			pass
		else:
			raise

	try:
		os.makedirs(os.path.join((app.config['PATH_TO_LSA']), prefix, suffix))
	except OSError:
		if os.path.isdir(os.path.join((app.config['PATH_TO_LSA']), prefix, suffix)):
			pass
		else:
			raise

	filename = str(prefix) + '/' + str(suffix) + '/' + "lsa_" + str(query) + "_" + str(k) + ".json"

	completeName = os.path.join(save_path, filename) #with the query for a name
	with open(completeName, 'w') as outfile:
		json.dump(jsonDict, outfile)

#Loads the LSA json for vis
#First checks if the file exists. If the file doesn't exist, it makes it.
def load_lsa(query, k):
	save_path = (app.config['PATH_TO_LSA'])
	pmid_list = query.split('+')  # list of string pmids
	pmid = pmid_list[0]  # get the first
	prefix = pmid[0:3]
	suffix = pmid[3:6]
	filename = str(prefix) + '/' + str(suffix) + '/' + "lsa_" + str(query) + "_" + str(k) + ".json"
	completeName = os.path.join(save_path, filename)
	try:
		with open(completeName) as infile:
			jsonDict = json.load(infile)
	except Exception as e:
		#there's no file! So we gotta make one!
		#load lemmas
		lemma_samples = load_lemma_cache(query)
		lemmas_for_lsa = [l[1] for l in lemma_samples]
		jsonDict = run_lsa1(lemmas_for_lsa, k)
		print_lsa(query, jsonDict, k)
	return jsonDict

###### L D A ########
def run_lda1(lda_lemmas, num_topics, n_top_words): #set at defulat k=3, number of words=5
	logging.info('Beginning Latent Dirichlet Allocation')
	tfidf, tfidf_vectorizer = get_tfidf(lda_lemmas)
	lda = fit_lda(tfidf, num_topics)
	jsonLDA = topics_lda(tfidf_vectorizer, lda, n_top_words)
	return jsonLDA

#Save the json for @app.route('/reslda/')
def print_lda(query, jsonLDA, k, w):
	save_path = (app.config['PATH_TO_LDA'])
	pmid_list = query.split('+')  # list of string pmids
	pmid = pmid_list[0]  # get the first
	prefix = pmid[0:3]
	suffix = pmid[3:6]

	try:
		os.makedirs(os.path.join((app.config['PATH_TO_LDA']), prefix))
	except OSError:
		if os.path.isdir(os.path.join((app.config['PATH_TO_LDA']), prefix)):
			pass
		else:
			raise

	try:
		os.makedirs(os.path.join((app.config['PATH_TO_LDA']), prefix, suffix))
	except OSError:
		if os.path.isdir(os.path.join((app.config['PATH_TO_LDA']), prefix, suffix)):
			pass
		else:
			raise

	filename = str(prefix) + '/' + str(suffix) + '/' + "lda_" + str(query) + "_" + str(k) + "_" + str(w) + ".json"

	completeName = os.path.join(save_path, filename)  # with the query for a name
	with open(completeName, 'w') as outfile:
		json.dump(jsonLDA, outfile)

def load_lda(query, k, w):
	save_path = (app.config['PATH_TO_LDA'])
	pmid_list = query.split('+')  # list of string pmids
	pmid = pmid_list[0]  # get the first
	prefix = pmid[0:3]
	suffix = pmid[3:6]
	filename = str(prefix) + '/' + str(suffix) + '/' + "lda_" + str(query) + "_" + str(k) + "_" + str(w) + ".json"
	completeName = os.path.join(save_path, filename)
	try:
		with open(completeName) as infile:
			jsonLDA = json.load(infile)
	except Exception as e:
		# there's no file! So we gotta make one!
		# load lemmas
		lemma_samples = load_lemma_cache(query)
		lemmas_for_lda = [l[1] for l in lemma_samples]
		jsonLDA = run_lda1(lemmas_for_lda, k, w)
		print_lda(query, jsonLDA, k, w)
	return jsonLDA

##### E M B E D D I N G  ######
#Input: query, top N desired bin, k clusters
#Output: prints csv for force directed graph
#NB: does not connect to DB!
def run_embeddings(query, k_clusters, top_n):
	logging.info("in run_embeddings function")#
	logging.info(query)
	words, tags = get_words_tags(query) #list of words/tags per doc
	transformed_sentence = transform_text(words, tags)
	npDict = chooseTopNPs(transformed_sentence)
	logging.info("done with npDict")
	logging.info(type(top_n))
	logging.info(type(k_clusters))
	if top_n == 100:
		logging.info("w=1-100")
		top = list(npDict.most_common(top_n))
	elif top_n == 200:
		logging.info("w=101-200")
		top_nps100 = list(npDict.most_common(100))
		top_nps200 = list(npDict.most_common(200))
		top = [item for item in top_nps200  if item not in top_nps100]
	elif top_n == 300:
		logging.info("w=201-300")
		top_nps200 = list(npDict.most_common(200))
		top_nps300 = list(npDict.most_common(300))
		top = [item for item in top_nps300  if item not in top_nps200]
	elif top_n == 400:
		logging.info("w=301-400")
		top_nps300 = list(npDict.most_common(300))
		top_nps400 = list(npDict.most_common(400))
		top = [item for item in top_nps400  if item not in top_nps300]
	else:
		top = list(npDict.most_common(top_n))
	logging.info("done with top NPs")
	logging.info("loaded model all the way!")
	matrix = getNPvecs(top, fasttext_model) #loaded as global variable
	logging.info("getting the matrix!")
	kmeans = KMeans(n_clusters=k_clusters, random_state=2).fit(matrix)
	results = list(zip(kmeans.labels_, top))
	embedding_json(results, query, k_clusters, top_n) #this saves it as a file
	logging.info("made json for embedding topic model")

def embedding_lookup(query, k_clusters, top_n):
	save_path = (app.config['PATH_TO_FGRAPHS'])
	pmid_list = query.split('+')  # list of string pmids
	pmid = pmid_list[0]  # get the first
	prefix = pmid[0:3]
	suffix = pmid[3:6]
	filename = str(prefix) + '/' + str(suffix) + '/' + 'fgraph_' + str(query) + '_' + str(k_clusters) + '_' + str(top_n) + '.json'
	completeName = os.path.join(save_path, filename)
	returnName = os.path.join('fgraphs', filename)
	if os.path.isfile(completeName): #if the file exists, all's well!
		logging.info("fgraph file exists! :)")
		pass
	else:
		logging.info("failed to load/find an existing fgraph. Make a new one.")
		run_embeddings(query, k_clusters, top_n)
	return returnName


