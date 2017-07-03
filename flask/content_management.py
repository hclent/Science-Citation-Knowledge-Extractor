from processors import * #pyProcessors
import os.path, time, re, logging, pickle, json, codecs, arrow
import operator
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



## Supporting functions for app.py

################## LOGGING #######################################################

#Create log
logging.basicConfig(filename='.app.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

############# PROCESSORS SERVER ##################################################

#Set a PROCESSORS_SERVER environment variable.
#It may take a minute or so to load the large model files.
def connect_to_Processors(port_num):
  path = '/home/hclent/anaconda3/envs/py34/lib/python3.4/site-packages/processors/processors-server.jar'
  api = ProcessorsAPI(port=port_num, jar_path=path, keep_alive=True, jvm_mem="-Xmx100G")
  logging.info('Connected to pyProcessors')
  #Initialize the bionlp annotator by passing it a doc
  init_doc = api.bionlp.annotate("The mitochondria is the powerhouse of the cell.")
  return api

################ WORD EMBEDDING MODEL ############## (Global var)

#This is really slow to load :/ Maybe shouldn't be a global var?
# Don't want to re-load for each analysis though
fasttext_model = load_model('/home/hclent/tmp/fastText/16kmodel.vec')

################### INPUT #########################################################

#User can enter in as many pubmed ids as they want into text box
#This method creates a list of them
def multiple_pmid_input(user_input):
	logging.info('cleaning user input')
	clean = re.sub('\,', ' ', user_input)
	ids = clean.split() #list of pmids
	return ids

################### DATABASE #####################################################
#If the pmid is NOT in the db, we also need to getSelfText and write that info to db
def scrape_and_write_Input(user_input):
	logging.info('retrieve Self text, and write to db')
	self_pmcid, self_abstract_check, self_article_check = getSelfText(user_input) #Entrez_IR function
	logging.info("PMICD: ")
	logging.info(self_pmcid)
	updateInputPapers(user_input, self_pmcid, self_abstract_check, self_article_check)  # put getSelfText into database

#input: user_input pmid
#output: list of dicts with annotation checks [{"pmcid": 123, "annotated": yes}]
def annotation_check(user_input):
	a_check = []
	pmc_ids = db_citation_pmc_ids(user_input) #Used to use getCitationIDs(user_input) here but updated to using my own db
	for citation in pmc_ids:
		#print(citation)
		annotationDict = {"pmcid": citation, "annotated": []}

		prefix = '/home/hclent/data/pmcids/' + str(citation[0:3])  # folder for first 3 digits of pmcid
		suffix = prefix + '/' + str(citation[3:6])  # folder for second 3 digits of pmcid nested in prefix
		filename = suffix + '/' + str(citation) + '.json'
		with open(filename) as data_file:
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
def new_or_copy_db(citation): #citation is a dict
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
def run_IR_not_db(user_input):
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
	scrape_and_write_Input(user_input)

	#Now retrieve citations
	logging.info("Get basic info about the citations")
	# Previously unseen pmcids only in allCitationsInfo.
	# Previously seen pmcids are copied to db for new pmid in getCitedInfo
	allCitationsInfo = getCitedInfo(pmc_ids, user_input) #output: list of dictionaries [{pmid: 1234, author: human, ...}]
	logging.info("Write basic citation info to citations db")
	for citation in allCitationsInfo:
		logging.info(citation)
		new_or_copy_db(citation)

	#Get content and update citations db
	contentDictList = getContentPMC(pmc_ids, user_input)
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
def run_IR_in_db(user_input):
	logging.info('PMID is in the database')
	# Check for new papers:
	num_in_db = db_input_citations_count(user_input) #checks MY db
	pmc_ids = getCitationIDs(user_input) #checks ENTREZ DB
	num_current = len(pmc_ids)
	#If there are new papers,
	if num_current == num_in_db: #TODO change this back to > after i've fixed authors problem
		need_to_annotate = 'yes'
		#print("there are new citations!", (num_current, num_in_db))
		#update number of citations in inputPaper db

		update = inputPapers.update().\
			where(inputPapers.c.pmid == user_input).\
			values(num_citations=num_current)
		conn.execute(update)

		#now get the new citation info
		allCitationsInfo = getCitedInfo(pmc_ids, user_input)  # output: list of dictionaries [{pmid: 1234, author: human, ...}] #skips duplicates
		logging.info("Write basic citation info to citations db for new papers")
		for citation in allCitationsInfo:
			new_or_copy_db(citation)

		#Get content and update citations db
		logging.info("now get the content for the new stuff")
		contentDictList = getContentPMC(pmc_ids, user_input)
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




#Depreciate this  :|
# def new_citations_from_db(user_input):
# 	apa_citations, db_journals, db_dates, db_urls = db_citations_retrieval(user_input)
# 	return apa_citations, db_urls
# 	#apa_citations called 'main' in app.py


#Data for populating statistics page in app
def get_statistics(query):
	total_pubs, unique_pubs, abstracts, whole, sentences, words =  db_query_statistics(query)
	statistics = [total_pubs, unique_pubs, abstracts, whole, sentences, words]
	return statistics



#use query to get info about input papers
def statsSelfInfo(query):
    input_click_citations = []
    pmid_list = query.split('+')  # list of string pmids
    for user_input in pmid_list:
        apa = db_inputPapers_retrieval(user_input)
        url = "https://www.ncbi.nlm.nih.gov/pubmed/"+str(user_input)
        href_label = (apa, url) #store apa and url as a tuple
        input_click_citations.append(href_label) #then append to list
    return input_click_citations


#take a query and generate x and y datapoints for pubs x year bar chart in "Stats" tab
#TODO: sepperate counts (have a different bar) for each inputPaper
def stats_barchart(query):
	pmid_list = query.split('+') #list of string pmids

	x_all = [] #list of list
	y_all = []
	for user_input in pmid_list:
		journals = []
		dates = []
		db_journals, db_dates = db_bar_chart(user_input)
		for j in db_journals:
			journals.append(j)
		for d in db_dates:
			dates.append(d)
		#TODO: only issue here is that we need the years range for the QUERY, and also 0 counts for any empties.
		#maybe something like:
		#years_range = getfromquerydb(qyer)
		#x, y = paper_dates_barchart(journals, dates, years_range, user_input)
		x, y = paper_dates_barchart(journals, dates, user_input)
		x_all.append(x)
		y_all.append(y)
	return x_all, y_all



############ DATA VISUALIZATIONS #################################################
#Updated to SqlAlchemy
#TODO: no mechanism for updating **DB** (table: queries) if more citations have been found!!
#TODO: update wc cache if new citations to a query
#TODO: check for journals vis file
def print_journalvis(query):
	record = checkForQuery(query)  # check for query in db.
	logging.info("checked for query!!!")
	logging.info(record)
	if record == 'empty':
		#if the record has never been seen before, do the journalsvis and write to db
		years_range = get_years_range(query) #need range for ALL journals, not just last one
		logging.info(years_range)
		publication_data, range_info = journals_vis(years_range, query) #range info = [('2008', '2016'), 165, 48]
		logging.info(range_info)
		journal_years = range_info[0]
		logging.info(journal_years)
		q = '+'
		logging.info(q)
		range_years = str(q.join(journal_years))
		logging.info("range years: "+range_years)
		unique_publications = range_info[1]
		unique_journals = range_info[2]
		logging.info(range_info)
		logging.info('Printing JOURNALS to JSON')
		save_path = '/home/hclent/data/journals/' #save in journals folder
		completeName = os.path.join(save_path, ('journals_'+(str(query))+'.json')) #named after query
		with open(completeName, 'w') as outfile:
			json.dump(publication_data, outfile)
		date = str(arrow.now().format('YYYY-MM-DD'))

		update = queries.insert().\
			values(dict(datestamp=date, query=query, range_years=range_years, unique_pubs=unique_publications,
						unique_journals=unique_journals))
		conn.execute(update)


	if record == 'yes': #if its in the db, just get the important things from the db!!
		range_years, unique_publications, unique_journals = getJournalsVis(query)
	return range_years, unique_publications, unique_journals



#TODO: chcek for wc file
#TODO: update wc cache if new citations to a query
def vis_wordcloud(neslist, nes_categories, w_number):
	nesDict = frequency_dict(neslist, nes_categories)
	wcl = wordcloud(nesDict, int(w_number))
	return wcl


def vis_heatmapTitles(lemma_samples, years):
	titles = []  # want citations instead of titles

	zipped = zip(years, lemma_samples)
	sorted_data = list(sorted(zipped, key=lambda z: int(z[0])))
	pmcids = [l[1][0] for l in sorted_data]

	#pmcids = [l[0] for l in lemma_samples]
	for id in pmcids:
		hyperlink = db_citations_hyperlink_retrieval(id)
		titles.append(hyperlink[0])
	#titles =
	return titles


def vis_heatmap(lemma_samples, nes_samples, nes_categories, w_number):
	neslist = [n[1] for n in nes_samples]
	nesDict = frequency_dict(neslist, nes_categories)
	#everything is sorted by year inside doHeatmap
	x_docs, y_words, z_counts, years = doHeatmap(nesDict, w_number, lemma_samples)
	#sorted by year in vis_HeatmapTitles the sme way as doHeatmap
	titles = vis_heatmapTitles(lemma_samples, years)
	return x_docs, y_words, z_counts, titles


#TODO: vis_heatmap and vis_clustermap bascially work the same.
#TODO: share the x_docs, y_words, z_counts between the two, for faster work :)
#TODO: fix papers axis so its not smooshed together and ugly
#TODO: fix the way papers are saved here, so we have the nested directories for faster lookup
#TODO: check for clustermap file
#TODO: update clustermap cache if new citations to a query
def vis_clustermap(lemma_samples, nes_samples, nes_categories, w_number, query):
	logging.info("starting clustermap")
	x, y, z, years = vis_heatmap(lemma_samples, nes_samples, nes_categories, w_number)
	logging.info("making clustermap data")
	seaData = make_seaborn_data(x, y, z)
	logging.info("saving clustermap png")
	saveName = makeClusterMap(seaData, query)
	return saveName #return filename


def vis_kmeans(lemma_samples, num_clusters):
	#use query to get titles
	titles = [] #want citations instead of titles

	pmcids = [l[0] for l in lemma_samples]
	for id in pmcids:
		hyperlink = db_citations_hyperlink_retrieval(id)
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
def inputEligible(query):
	papers = []
	values = ['paper1', 'paper2', 'paper3', 'paper4', 'paper5']
	path_to_paper = []
	pmid_list = query.split('+')  # list of string pmids
	display_title = [] #title for the dropdown menu

	for pmid in pmid_list:
	#We are assuming that if there is a pmcid for the pmid, that means it was in PubmedCentral and we scraped it
		pmcid = pmid2pmcid(pmid)
		if pmcid != "NA":
			#print(pmcid)
			#print(pmid + " = " + pmcid)
			#get the pmcid of the pmid

			prefix = pmcid[0:3]
			suffix = pmcid[3:6]
			filename = '/home/hclent/data/pmcids/' + str(prefix) + '/' + str(suffix) + '/' + str(pmcid) + '.txt' # look in folder that matches pmcid
			#print(filename)
			truth_value = os.path.isfile(filename)
			if truth_value is True:
				#print(filename) #NA won't exist
				papers.append(pmid)
				path_to_paper.append(filename)
				display = db_pmid_axis_label(pmid)
				keep_display = display[0]
				display_title.append(keep_display)

	eligible_papers = list(zip(values, papers, path_to_paper, display_title))
	return eligible_papers


#visualization for scifi div
def vis_scifi(corpus, query, eligible_papers):
	corpus_vec, color = load_corpus(corpus, eligible_papers)
	eligible_cosines = get_cosine_eligible(corpus_vec, eligible_papers)
	logging.info("eligible_cosines: done")
	data_vecs_list, pmcids_list  = load_datasamples(query)
	logging.info("data_vecs_list: done")
	cosine_list = get_cosine_list(corpus_vec, data_vecs_list)
	logging.info("cosine_list: done")
	#cosine list is the cosine sim score for each document
	sorted_combos = add_urls(cosine_list, color, pmcids_list)
	logging.info("sorted combos: done")
	all_sorted_combos = add_eligible_cosines(sorted_combos, eligible_papers, eligible_cosines)
	logging.info("all_sorted_combos: done")
	x, y, names, color_list = prepare_for_histogram(all_sorted_combos)
	logging.info(x)
	logging.info(y)
	logging.info(color_list)
	return x, y, names, color_list

# eligible_papers = [('paper1', '18952863', '/home/hclent/data/pmcids/259/367/2593677.txt')]
# corpus = 'darwin'
# query = "18952863+18269575"
# x, y, names, color_list = vis_scifi(corpus, query, eligible_papers)




############ PROCESSING BIODOCS ############################################
#Take pmcid.txt and get an annotated document, as well as lemmas and named entities
#Doesn't re-annotated documents that have already been annotated.
#Updated to sqlalchemy
def do_multi_preprocessing(user_input):
	logging.info('Beginning multiprocessing for NEW (unprocessed) docs')
	t1 = time.time()
	docs = retrieveDocs(user_input)
	multiprocess(docs) #if docs is empty [], this function just passes :)
	# # #Now update annotated_check
	a_check = annotation_check(user_input)

	#TODO: only update/annotate check if there's new docs since last time?
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
	biodocs = retrieveBioDocs(user_input)
	biodoc_data = loadBioDoc(biodocs) #list of dictionaries[{pmid, lemmas, nes, sent_count, token_count}]
	#No problem getting biodocs or biodoc_data ... problem comes with updating db...
	#update db with sents and tokens
	for b in biodoc_data:
		update_annotations(b, user_input)
	logging.info("Execute everything: done in %0.3fs." % (time.time() - t1))
	return biodoc_data



############ TOPIC MODELING ############################################
def run_lsa1(lsa_lemmas, k):
	logging.info('Beginning Latent Semantic Analysis')
	tfidf, tfidf_vectorizer = get_tfidf(lsa_lemmas)
	jsonDict = do_LSA(tfidf, tfidf_vectorizer, k) #need to make this an option
	return jsonDict


def run_lda1(lda_lemmas, num_topics, n_top_words): #set at defulat k=3, number of words=5
	logging.info('Beginning Latent Dirichlet Allocation')
	tfidf, tfidf_vectorizer = get_tfidf(lda_lemmas)
	lda = fit_lda(tfidf, num_topics)
	jsonLDA = topics_lda(tfidf_vectorizer, lda, n_top_words)
	return jsonLDA


#Input: query, top N desired bin, k clusters
#Output: prints csv for force directed graph
def run_embeddings(query, top_n, k_clusters):
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
	embedding_json(results, query)
	logging.info("made json for embedding topic model")



########### WRITING TO JSON / PICKLE ###############################################

def print_lsa(query, user_input, jsonDict):
	#Save the json for @app.route('/reslsa/')
	logging.info('Printing LSA to JSON')
	save_path = '/home/hclent/data/topics/lsa' #in the folder of the last pmid
	completeName = os.path.join(save_path, ('lsa_'+(str(query))+'.json')) #with the query for a name
	with open(completeName, 'w') as outfile:
		json.dump(jsonDict, outfile)

def print_lda(query, user_input, jsonLDA):
	#Save the json for @app.route('/reslda/')
	logging.info('Printing LDA to JSON')
	save_path = '/home/hclent/data/topics/lda' #in the folder of the last pmid
	completeName = os.path.join(save_path, ('lda_'+(str(query))+'.json'))  #with the query for a name
	with open(completeName, 'w') as outfile:
		json.dump(jsonLDA, outfile)


def flatten(listOfLists):
    return list(chain.from_iterable(listOfLists))


############## GRAVEYARD ##########################################################
