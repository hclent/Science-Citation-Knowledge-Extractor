from flask import Flask, render_template, request, flash, url_for, redirect, session, g, Blueprint
from flask_wtf import Form
from wtforms import TextField, SelectField
import gc, time, datetime, pickle, os.path
import sys, csv
from werkzeug.serving import run_simple
sys.path.append('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/')
from database_management import engine, connection, inputPapers #mine
from content_management import * #mine
from citation_venn import make_venn #mine
from processors import *


app = Flask(__name__, static_url_path='/hclent/Webdev-for-bioNLP-lit-tool/flask/static')
app.config.from_pyfile('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/configscke.cfg', silent=False) #pass abs path


#Create Form for handling user-entered pmid
#Need to pass pmid in form to Entrez_IR.py
class pmidForm(Form):
	pmid = TextField('PubmedID')


#Main page
#Prints sample results from 2 coge publications
#User inputs a pubmed id and is then redirected to /results
@app.route("/cogecrawl/")
def cogecrawl():
	#error = None
	with open('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_citations.pickle', 'rb')as f:
		citations_with_links = pickle.load(f)
	return render_template("dashboard.html", citations_with_links=citations_with_links)



#Getting Flask-WTFs to work with sqlite3 here
#This function uses user entry to run the Entrez_IR.py 
#User entered pmid is entered into sqlite3 database
@app.route('/results/', methods=["POST"])
def results():
	logging.info("In app route RESULTS")
	form = pmidForm()
	try:
		if request.method == 'POST':
			entry = form.pmid.data #THIS IS THE USER INPUT FROM THE FORM #referencing 'class pmidForm'
			pmid_list = multiple_pmid_input(entry) #list for handling multiple pmids
			logging.info(pmid_list)


			# If the user inputs more than 5 PMIDs, return the home page and flash a warning
			# Need 5 PMIDs or less
			if len(pmid_list) > 5:
				flash('You have entered more than 5 PMIDs. Please reduce your query to 5 PMIDs or less to continue.')
				with open('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_citations.pickle', 'rb')as f:
					citations_with_links = pickle.load(f)
				return render_template("dashboard.html", citations_with_links=citations_with_links)


			q = '+'
			query = str(q.join(pmid_list))
			logging.info("query: " + str(query))

			main_info = [] #main_info and target_urls for "citaitons" page
			target_urls = []

			#TODO: need to re-factor how i do data_samples and ners as DICTS with pmcid keys
			data_samples = []
			ners = []


			for user_input in pmid_list:
				logging.info(str(user_input))
				user_input = str(user_input)

				############################################
				#Check database for pmid #Does the entry exists in the db already?
				conn = connection()
				s = inputPapers.select().\
					where(inputPapers.c.pmid == user_input)
				c = conn.execute(s)
				check1 = c.fetchone()
				c.close()


				#if the entry does NOT exist in the db already, will need to retrieve text and annotate
				if check1 is None:
					flash('new pubmedid!')
					#Using user_input for Information Retireval of citing pmcids and info about them
					run_IR_not_db(user_input)
					logging.info("beginning multi-preprocessing")
					biodoc_data = do_multi_preprocessing(user_input)
					logging.info("done with new document multi_preprocessing")
					logging.info("writing the BIODOC LEMMAS")
					biodoc_to_db(biodoc_data) #writes data_samples (lemmas) and NER to db
					logging.info("* done writing the biodoc lemmas")
					#TODO: Change how I am making data_samples. Use dict so I can access by pmcid, NOT just index
					#TODO: Change how I am making named_entites. Use dict so I can access by pmcid, NOT just index


					#After all citations have been processed, now we can do the analyses:
					if user_input == pmid_list[-1]: #if its the last pmid
						logging.info("last pmid in the query")
						logging.info("begin journal vis!")


						#JOURNALS VIS STUFF HERE
						logging.info(user_input+" is the last one (JOURNALS)")
						range_years, unique_publications, unique_journals = print_journalvis(query)#TODO: Check for vis.json first

						#TOPIC MODELING HERE
						# logging.info(user_input+" is the last one (LSA)")
						# #Do Latent Semantic Analysis and return jsonDict for data vis
						# jsonDict = run_lsa1(data_samples, 2)
						# print_lsa(query, user_input, jsonDict) #print lsa topic model to json
                        #
						# logging.info(user_input+" is the last one (LDA)")
						# jsonLDA = run_lda1(data_samples, 3, 5)
						# print_lda(query, user_input, jsonLDA) #print lda topic model to json


						## FUNCTION THAT CACHES DATA_SAMPLES AND NES HERE:
						#print_data_and_nes(query, user_input, data_samples, ners) #print data_samples and nes_list to pickle


					#after info written to db, now can access db and get formated main_info (main)
					#TODO: why is this here? and not somewhere else in the code?
					main, db_urls = new_citations_from_db(user_input)
					for mi in main:
						main_info.append(mi)
					logging.info("done with main info list")
					for url in db_urls:
						target_urls.append(url)
					logging.info("done with url list")




				#if the entry IS in the db, no need to retrieve text from Entrez, just grab from db
				if check1 is not None:
					flash("alreay exists in database :) ")
					#Using user_input for Information Retireval - checks if any new papers have been added that we need to scrape
					need_to_annotate = run_IR_in_db(user_input)
					biodoc_data = do_multi_preprocessing(user_input)
					logging.info("done with new document multi_preprocessing")
					logging.info("writing the BIODOC LEMMAS")
					#biodoc_to_db(biodoc_data)  # writes data_samples (lemmas) and NER to db ## writing to db suucckkss
					print_data_samples(user_input, biodoc_data)
					# if need_to_annotate == 'yes':
					# 	logging.info("need to annotate new documents")
					# 	biodoc_data = do_multi_preprocessing(user_input)
					# 	biodoc_to_db(biodoc_data)
					#   print_data_samples(user_input, biodoc_data)
					# if need_to_annotate == 'no':
					# 	logging.info("dont need to annotate any new documents")
					# 	pass

					main, db_urls = new_citations_from_db(user_input)
					for mi in main:
						main_info.append(mi)
					logging.info("done with main info list")
					for url in db_urls:
						target_urls.append(url)
					logging.info("done with url list")

					## Now that we have all the data, do the topic model
					## Only want to save final topic model (not running topic model)
					if user_input == pmid_list[-1]:
						logging.info("last pmid in the query")
						logging.info("begin journal vis!")

						# JOURNALS VIS STUFF HERE
						logging.info(user_input + " is the last one (JOURNALS)")
						range_years, unique_publications, unique_journals = print_journalvis(query)  # TODO: Check for vis.json first


						## PRINT Data_samples and stuff
						#TODO: Depreciated
						#print_data_and_nes(query, user_input, data_samples, nes_list)

						### TOPIC MODELING STUFF. SHOULD DO ALL IN iFRAMES
						# logging.info(user_input + " is the last one (LSA)")
						# jsonDict = run_lsa1(data_samples, 2)
						# print_lsa(query, user_input, jsonDict)  # print lsa topic model to json
                        #
						# logging.info(user_input + " is the last one (LDA)")
						# jsonLDA = run_lda1(data_samples, 3, 5)
						# print_lda(query, user_input, jsonLDA)  # print


				#Housekeeping
				gc.collect() #garbage collector for cleaning up unneeded stuff
				session['entered_id'] = True
				session['engaged'] = 'engaged'



		citations_with_links = list(zip(main_info, target_urls))


		return render_template('results.html', form=form, citations_with_links=citations_with_links,
	   			main_info = main_info,  query=query, range_years=range_years, unique_publications=unique_publications, unique_journals=unique_journals)
				#took out start_year, end_year

	except Exception as e:
		return(str(e))




################ Forms for visualization toggle ################
#Form for Topic Models
class visOptions(Form):
	k_val = SelectField('k_val', choices=[(2,'k=2'),(3,'k=3'),(4,'k=4'),(5,'k=5'),(6, 'k=6'),
										  (7, 'k=7'),(8, 'k=8'),(9, 'k=9'),(10,'k=10'),
										  (11, 'k=11'),(12, 'k=12'),(13, 'k=13')])
	w_words = SelectField('w_words', choices=[(4, 'w=4'),(5, 'w=5'),(6, 'w=6'), (7, 'w=7'),
											  (100, 'w=1-100'),(200, 'w=101-200'),(300, 'w=201-300'),(400, 'w=301-400')])


#Form for NEs tab
class nesOptions(Form):
	w_words = SelectField('w_words', choices=[(2, 'N'),(3, '3'),(10, '10'), (25, '25'), (50, '50'),(100, '100'),(200, '200'),
											  (300, '300')])

#Form for TextCompare tab
class corpusOptions(Form):
	corpus = SelectField('corpus', choices=[('startrek', 'startrek'),('frankenstein', 'frankenstein'),('youth', 'youth'),
											('darwin', 'darwin'), ('austen', 'austen'),
											('brain_speech','brain_speech'), ('bible','bible'), ('grecoroman','grecoroman'),
											('last_evolution','last_evolution'), ('mars','mars'), ('mouse','mouse'),
											('sherlock','sherlock'),('yeast','yeast'),
											('paper1', 'paper1'),
											('paper2', 'paper2'),('paper3', 'paper3'), ('paper4', 'paper4'), ('paper5', 'paper5')])

################ Default CoGe Data #############################
@app.route('/cogembeddings/', methods=["GET","POST"]) #default coge embeddings topic for iframe
def cogeembeddings():
	form = visOptions()
	if request.method =='POST':
		logging.info("posted something to cogembeddings")
		query = '18952863+18269575'
		window = int(form.w_words.data)
		logging.info(window)
		k_clusters = int(form.k_val.data)  # 2,3,4,or 5
		logging.info(k_clusters)
		run_embeddings(query, window, k_clusters)  # 50 words in 6 clusters
		filepath = 'fgraphs/fgraph_' + str(query) + '.json'
		logging.info(filepath)
		return render_template('coge_embeddings.html', filepath=filepath)
	else:
		filepath = 'coge_embed.json'
		return render_template('coge_embeddings.html', filepath=filepath)

@app.route('/cogelsa/', methods=["GET","POST"]) #default coge lsa for iframe
def cogelsa():
	form = visOptions()
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		query = '18952863+18269575'

		#UPDATING TO NEW & IMPROVED lemma_samples!
		lemma_file = '/home/hclent/data/pmcids/189/528/lemma_samples_18952863+18269575.pickle'
		with open(lemma_file, "rb") as f:
			lemma_samples = pickle.load(f)

		lemmas_for_lsa = [l[1] for l in lemma_samples] #ignore the pmcid's in l[0]

		logging.info("rerunning the analysis")
		k = int(k_clusters)

		jsonLSA = run_lsa1(lemmas_for_lsa, k)
		logging.info("did it all!")
		return render_template('coge_lsa.html', form=form, jsonLSA=jsonLSA) #needs to be parsed
	else: #if nothing is
		completeName = "/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_lsa.json"
		with open(completeName) as load_data:
			jsonLSA = json.load(load_data) #doesn't need to be parsed but unsure how to write that in javascript
		#so i'm going to read it in as a string that needs to be parsed anyway
		jsonLSA = re.sub('\'', '\"', str(jsonLSA)) #json needs double quotes, not single quotes
		return render_template('coge_lsa.html', form=form, jsonLSA=jsonLSA)


@app.route('/cogelda/', methods=["GET","POST"]) #default coge lda for iframe
def cogelda():
	form = visOptions()
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		num_words = form.w_words.data
		logging.info("the w value is "+str(num_words))
		query = '18952863+18269575'

		lemma_file = '/home/hclent/data/pmcids/189/528/lemma_samples_18952863+18269575.pickle'
		with open(lemma_file, "rb") as f:
			lemma_samples = pickle.load(f)

		lemmas_for_lda = [l[1] for l in lemma_samples]  # ignore the pmcid's in l[0]

		logging.info("rerunning the analysis")
		k = int(k_clusters)
		w = int(num_words)
		jsonLDA = run_lda1(lemmas_for_lda, k, w)
		return render_template('coge_lda.html', form=form, jsonLDA=jsonLDA)
	else:
		completeName = "/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_lda1.json"
		with open(completeName) as load_data:
			jsonLDA = json.load(load_data) #doesn't need to be parsed but unsure how to write that in javascript
		#so i'm going to read it in as a string that needs to be parsed anyway
		jsonLDA = re.sub('\'', '\"', str(jsonLDA)) #json needs double quotes, not single quotes
		return render_template('coge_lda.html', form=form, jsonLDA=jsonLDA)


@app.route('/cogejournals/') #default coge journals for iframe
def cogejournals():
	with open("/home/hclent/data/journals/journals_18952863+18269575.json") as load_data:
		journals = json.load(load_data)
	return render_template('coge_journals.html', journals=journals)


@app.route('/cogewordcloud/', methods=["GET","POST"]) #default coge NES Word Cloud for iframe
def cogewordcloud():
	form = nesOptions()
	if request.method == 'POST':

		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		nes_list =  pickle.load(open("/home/hclent/data/nes/nes_18952863+18269575.pickle", "rb")) #pre-processed already

		#logging.info(nes_list)
		wordcloud_data = vis_wordcloud(nes_list, nes_categories, w_number)
		popup = ' '
		return render_template('coge_wordcloud.html',  wordcloud_data=wordcloud_data, popup=popup)
	else:
		#Default data
		completeName = "/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_wcloud1.json"
		with open(completeName) as load_data:
			wordcloud_data = json.load(load_data) #this result doesn't need to be parsed but unsure how to write that in javascript
		#so i'm going to read it in as a string that needs to be parsed anyway
		wordcloud_data = re.sub('\'', '\"', str(wordcloud_data)) #json needs double quotes, not single quotes
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong> Default: N=10, from all categories.</div>'
		return render_template('coge_wordcloud.html', wordcloud_data=wordcloud_data, popup=popup)


@app.route('/cogeheatmap/', methods=["GET","POST"]) #default coge NES heatmap for iframe
def cogeheatmap():
	form = nesOptions()
	query = '18952863+18269575'
	if request.method == 'POST':
		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		nes_list =  pickle.load(open("/home/hclent/data/nes/nes_18952863+18269575.pickle", "rb")) #pre-processed already
		data_samples =  pickle.load(open("/home/hclent/data/data_samples/data_samples_18952863+18269575.pickle", "rb")) #pre-processed already

		x_docs, y_words, z_counts = vis_heatmap(data_samples, nes_list, nes_categories, w_number)
		titles = vis_heatmapTitles(query)
		#print(z_counts)
		#print(x_docs)
		#print(y_words)
		return render_template('coge_heatmap2.html', z_counts=z_counts, x_docs=x_docs, y_words=y_words, titles=titles)
	else:
		#Default data
		return render_template('coge_heatmap1.html')


@app.route('/cogeclustermap/', methods=["GET","POST"]) #default coge clustermap
def cogeclustermap():
	query = '18952863+18269575'
	form = nesOptions()
	if request.method == 'POST':
		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is " + str(w_number))
		nes_list = pickle.load(open("/home/hclent/data/nes/nes_18952863+18269575.pickle", "rb"))  # pre-processed already
		data_samples = pickle.load(open("/home/hclent/data/data_samples/data_samples_18952863+18269575.pickle", "rb"))  # p
		saveName = vis_clustermap(data_samples, nes_list, nes_categories, w_number, query)
		image = '/images/'+saveName
		return render_template('coge_clustermap.html', image=image)
	else:
		image = "images/coge_clustermap.png"
		return render_template('coge_clustermap.html', image=image)


@app.route('/cogekmeans/', methods=["GET","POST"]) #default coge k-means clustering for iframe
def cogekmeans():
	form = visOptions()
	if request.method == 'POST':

		#depreciate data_samples
		#data_samples =  pickle.load(open("/home/hclent/data/data_samples/data_samples_18952863+18269575.pickle", "rb")) #pre-processed already
		lemma_file = '/home/hclent/data/pmcids/189/528/lemma_samples_18952863+18269575.pickle'
		with open(lemma_file, "rb") as f:
			lemma_samples = pickle.load(f) #the pmcids are in lemma_samples yay!

		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates, titles0, titles1, titles2, titles3, titles4 = vis_kmeans(lemma_samples, k_clusters)

		return render_template('coge_kmeans2.html', x0_coordinates=x0_coordinates, y0_coordinates=y0_coordinates, z0_coordinates=z0_coordinates,
							   x1_coordinates=x1_coordinates, y1_coordinates=y1_coordinates, z1_coordinates=z1_coordinates,
							   x2_coordinates=x2_coordinates, y2_coordinates=y2_coordinates, z2_coordinates=z2_coordinates,
							   x3_coordinates=x3_coordinates, y3_coordinates=y3_coordinates, z3_coordinates=z3_coordinates,
							   x4_coordinates=x4_coordinates, y4_coordinates=y4_coordinates, z4_coordinates=z4_coordinates,
							   titles0=titles0, titles1=titles1, titles2=titles2, titles3=titles3, titles4=titles4)

	else:
		return render_template('coge_kmeans.html')



@app.route('/coge_stats/') #default coge statistics for iframe
def coge_stats():
	query = "18952863+18269575"
	input_click_citations = statsSelfInfo(query)
	return render_template('coge_stats.html', input_click_citations=input_click_citations)


@app.route('/coge_scifi/', methods=["GET","POST"]) #default coge scifi for iframe
def coge_scifi():
	form = corpusOptions()
	eligible_papers = [('paper1', '18952863', '/home/hclent/data/pmcids/259/367/2593677.txt')]
	if request.method == 'POST':
		logging.info("posted a thing in scifi!")
		corpus = form.corpus.data
		logging.info(corpus)
		query = "18952863+18269575"
		if corpus == 'darwin':
			title = 'On The Origin of Species'
		if corpus == 'yeast':
			title = 'Yeast by Thomas Henry Huxley'
		if corpus == 'mouse':
			title = 'The Dancing Mouce, a Study in Anerimal Behavior by Robert Yerkes'
		if corpus == 'brain_speech':
			title = 'The Brain and The Voice in Speech and Song by Mott Frederick Walker'
		if corpus == 'grecoroman':
			title = 'Outlines of Greek and Roman Medicine by James Sands Elliott'
		if corpus == 'startrek':
			title = 'Star Trek: The Next Generation'
		if corpus == 'mars':
			title = 'Guilliver of Mars by Edwin Arnold'
		if corpus == 'last_evolution':
			title = 'The Last Evolution by John W. Campbell'
		if corpus == 'youth':
			title = 'Youth by Isaac Asimov'
		if corpus == 'frankenstein':
			title = 'Frankenstein; or, The Modern Prometheus'
		if corpus == 'sherlock':
			title = 'Sherlock Holmes by Sir Arthur Conan Doyle'
		if corpus == 'austen':
			title = 'Pride and Prejudice'
		if corpus == 'bible':
			title = 'The Bible'
		if corpus == 'paper1':
			title =  'PMID: 18952863'
		x, y, names, color = vis_scifi(corpus, query, eligible_papers)
		return render_template('coge_scifi2.html', x=x, y=y, title=title, color=color, names=names, eligible_papers=eligible_papers)
	else:
		flash('Some input paper(s) are not avaliable for TextCompare')
		return render_template('coge_scifi.html', eligible_papers=eligible_papers)


############### Results visualizations #########################################
@app.route('/resjournals/<query>/<range_years>', methods=["GET", "POST"]) #user journals for iframe
def resjournals(query, range_years):
	#need to get last user_input
	logging.info("in routine res-journals")
	logging.info("JOURNALS ID: " +str(query))
	logging.info("YEARS RANGE: " +str(range_years))
	#Need years for range
	years_list = range_years.split('+')
	s_year = years_list[0]
	e_year = years_list[1]

	#only want to load the json for the LAST id in the query (so includes all)
	pmid_list = query.split('+') #list of string pmids
	last_entry = pmid_list[-1]
	logging.info("the last entry is: " + str(last_entry))

	file_name = "journals_"+str(query)+".json"
	logging.info("last entry's JOURNAL is named: " + str(file_name))
	savePath = "/home/hclent/data/journals"
	completeName = os.path.join(savePath, file_name)
	logging.info("complete file: " + str(completeName))
	with open(completeName) as load_data:
		journals = json.load(load_data)
	return render_template('results_journals.html', journals=journals, s_year=s_year, e_year=e_year)



#TODO: re-implement resembeddings for results!
@app.route('/resembed/<query>', methods=["GET", "POST"]) #user embeddings for iframe
def resembeddings(query):
	if request.method == 'POST':
		return render_template('results_embeddings.html')
	else:
		return render_template('results_embeddings.html')



@app.route('/reslsa/<query>', methods=["GET", "POST"]) #user lsa for iframe
def reslsa(query):
	form = visOptions()
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))

		data_filename = "/home/hclent/data/data_samples/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		num_pubs = int(len(data_samples))
		logging.info("there are  "+str(num_pubs)+ " publications")
		logging.info("rerunning the analysis")
		k = int(k_clusters)
		if num_pubs < k:
			logging.info("k value is larger than number of publications")
			#flash("For LSA, you cannot have more topics than documents. Try again")
		temp_jsonDict = run_lsa1(data_samples, k)
		logging.info("did it all!")
		return render_template('results_lsa.html', query=query, jsonDict=temp_jsonDict)
	else:
		#need to get last user_input
		#use id to do stuff
		logging.info("in routine resla")
		logging.info("RES-LSA ID: " +str(query))
		#only want to load the json for the LAST id in the query (so includes all)
		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		logging.info("the last entry is: " + str(last_entry))
		file_name = "lsa_"+str(query)+".json" #file named after query
		logging.info("last entry's LSA is named: " + str(file_name))
		savePath = "/home/hclent/data/topics/lsa" #saved in folder of topics/lsa
		completeName = os.path.join(savePath, file_name)
		logging.info("complete file: " + str(completeName))
		with open(completeName) as load_data:
			jsonDict = json.load(load_data)
		return render_template('results_lsa.html', query=query, jsonDict=jsonDict)


@app.route('/reslda/<query>', methods=["GET", "POST"]) #user lda for iframe
def reslda(query):
	form = visOptions()
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		num_words = form.w_words.data
		logging.info("the w value is "+str(num_words))

		data_filename = "/home/hclent/data/data_samples/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		logging.info("rerunning the analysis")
		k = int(k_clusters)
		w = int(num_words)
		temp_jsonLDA = run_lda1(data_samples, k, w)
		return render_template('results_lda.html', form=form, jsonLDA=temp_jsonLDA, query=query)
	else:
		#need to get last user_input
		#use id to do stuff
		logging.info("in routine reslDa")
		logging.info("RES-LDA ID: " +str(query))
		file_name = "lda_"+str(query)+".json" #file named after query
		logging.info("last entry's LDA is named: " + str(file_name))
		savePath = "/home/hclent/data/topics/lda" #saved in folder of topics/lda
		completeName = os.path.join(savePath, file_name)
		#logging.info("complete file: " + str(completeName))
		with open(completeName) as load_data:
			jsonLDA = json.load(load_data)
		logging.info(jsonLDA)
		return render_template('results_lda.html', form=form, jsonLDA=jsonLDA, query=query)


@app.route('/reswordcloud/<query>', methods=["GET", "POST"]) #user wordcloud for iframe
def reswordcloud(query):
	form = nesOptions()
	if request.method == 'POST':

		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))


		filename = "/home/hclent/data/nes/nes_"+str(query)+".pickle"
		nes_list =  pickle.load(open(filename, "rb")) #pre-processed already

		#print(nes_list)
		wordcloud_data = vis_wordcloud(nes_list, nes_categories, w_number)
		popup = ' '
		return render_template('results_wordcloud.html', query=query,  wordcloud_data=wordcloud_data, popup=popup)
	else:
		nes_categories= ['BioProcess', 'CellLine', 'Cellular_component', 'Family', 'Gene_or_gene_product', 'Organ', 'Simple_chemical', 'Site', 'Species', 'TissueType']
		logging.info(nes_categories)
		w_number = 10
		logging.info("the w value is "+str(w_number))


		filename = "/home/hclent/data/nes/nes_"+str(query)+".pickle"
		nes_list =  pickle.load(open(filename, "rb")) #pre-processed already

		wordcloud_data = vis_wordcloud(nes_list, nes_categories, w_number)
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong> Default: N=10, from all categories.</div>'

		return render_template('results_wordcloud.html', query=query, wordcloud_data=wordcloud_data, popup=popup)


@app.route('/res_heatmap/<query>', methods=["GET", "POST"]) #user heatmap for iframe
def res_heatmap(query):
	form = nesOptions()
	if request.method == 'POST':
		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		nes_filename = "/home/hclent/data/nes/nes_"+str(query)+".pickle"
		nes_list =  pickle.load(open(nes_filename, "rb")) #pre-processed already

		data_filename = "/home/hclent/data/data_samples/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		x_docs, y_words, z_counts = vis_heatmap(data_samples, nes_list, nes_categories, w_number)
		titles = vis_heatmapTitles(query)
		# print(z_counts)
		# print(x_docs)
		# print(y_words)
		popup = ' '
		return render_template('results_heatmap.html', query=query, z_counts=z_counts, x_docs=x_docs, y_words=y_words, popup=popup, titles=titles)
	else:
		nes_categories= ['BioProcess', 'CellLine', 'Cellular_component', 'Family', 'Gene_or_gene_product', 'Organ', 'Simple_chemical', 'Site', 'Species', 'TissueType']
		#print(nes_categories)
		w_number = 10
		logging.info("the w value is "+str(w_number))


		nes_filename = "/home/hclent/data/nes/nes_"+str(query)+".pickle"
		nes_list =  pickle.load(open(nes_filename, "rb")) #pre-processed already

		data_filename = "/home/hclent/data/data_samples/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		x_docs, y_words, z_counts = vis_heatmap(data_samples, nes_list, nes_categories, w_number)
		titles = vis_heatmapTitles(query)
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong> Default: N=10, from all categories.</div>'
		return render_template('results_heatmap.html', query=query, z_counts=z_counts, x_docs=x_docs, y_words=y_words, popup=popup, titles=titles)



@app.route('/res_clustermap/<query>', methods=["GET", "POST"]) #user heatmap for iframe
def res_clustermap(query):
	form = nesOptions()
	if request.method == 'POST':
		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		nes_filename = "/home/hclent/data/nes/nes_"+str(query)+".pickle"
		nes_list =  pickle.load(open(nes_filename, "rb")) #pre-processed already

		data_filename = "/home/hclent/data/data_samples/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		saveName = vis_clustermap(data_samples, nes_list, nes_categories, w_number, query)
		image = '/images/' + saveName

		return render_template('results_clustermap.html',image=image, query=query)
	else:
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong> Choose N and categories to run clustermap.</div>'
		return render_template('results_clustermapH.html', query=query, popup=popup)


@app.route('/res_kmeans/<query>', methods=["GET", "POST"]) #user k-means for iframe
def res_kmeans(query):
	form = visOptions()
	if request.method == 'POST':

		pmid_list = query.split('+') #list of string pmids
		data_filename = "/home/hclent/data/data_samples/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates, titles0, titles1, titles2, titles3, titles4 = vis_kmeans(data_samples, k_clusters, pmid_list)
		return render_template('res_kmeans1.html', query=query,
		   x0_coordinates=x0_coordinates, y0_coordinates=y0_coordinates, z0_coordinates=z0_coordinates,
		   x1_coordinates=x1_coordinates, y1_coordinates=y1_coordinates, z1_coordinates=z1_coordinates,
		   x2_coordinates=x2_coordinates, y2_coordinates=y2_coordinates, z2_coordinates=z2_coordinates,
		   x3_coordinates=x3_coordinates, y3_coordinates=y3_coordinates, z3_coordinates=z3_coordinates,
		   x4_coordinates=x4_coordinates, y4_coordinates=y4_coordinates, z4_coordinates=z4_coordinates,
			titles0=titles0, titles1=titles1, titles2=titles2, titles3=titles3, titles4=titles4)
	else:
		return render_template('res_kmeans1.html', query=query)


@app.route('/res_stats/<query>', methods=["GET"]) #user statistics for iframe
def res_stats(query):
	pmid_list = query.split('+') #list of string pmids
	input_click_citations = statsSelfInfo(query)
	venn_data = make_venn(pmid_list)
	statistics = get_statistics(pmid_list)
	sum_total = statistics[0]
	unique = statistics[1]
	sum_abstracts = statistics[2]
	sum_whole = statistics[3]
	sum_sents = statistics[4]
	sum_tokens = statistics[5]
	#get x, y coordinates for pubs x year bar chart.
	x, y = stats_barchart(query)
	return render_template('results_stats.html', input_click_citations=input_click_citations,
						   venn_data=venn_data, sum_total=sum_total,
						   unique=unique, sum_abstracts=sum_abstracts, sum_whole=sum_whole,
						   sum_sents=sum_sents, sum_tokens=sum_tokens, x=x, y=y)


@app.route('/results_scifi/<query>', methods=["GET","POST"]) #default coge scifi for iframe
def results_scifi(query):
	form = corpusOptions()
	pmid_list = query.split('+')  # list of string pmids
	#decide eligible papers:
	eligible_papers = inputEligible(query)
	logging.info("eligible papers: " +str(eligible_papers))
	if request.method == 'POST':
		logging.info("posted a thing in scifi!")
		corpus = form.corpus.data
		logging.info(corpus)
		if corpus == 'darwin':
			title = 'On The Origin of Species'
		if corpus == 'yeast':
			title = 'Yeast by Thomas Henry Huxley'
		if corpus == 'mouse':
			title = 'The Dancing Mouce, a Study in Anerimal Behavior by Robert Yerkes'
		if corpus == 'brain_speech':
			title = 'The Brain and The Voice in Speech and Song by Mott Frederick Walker'
		if corpus == 'grecoroman':
			title = 'Outlines of Greek and Roman Medicine by James Sands Elliott'
		if corpus == 'startrek':
			title = 'Star Trek: The Next Generation'
		if corpus == 'mars':
			title = 'Guilliver of Mars by Edwin Arnold'
		if corpus == 'last_evolution':
			title = 'The Last Evolution by John W. Campbell'
		if corpus == 'youth':
			title = 'Youth by Isaac Asimov'
		if corpus == 'frankenstein':
			title = 'Frankenstein; or, The Modern Prometheus'
		if corpus == 'sherlock':
			title = 'Sherlock Holmes by Sir Arthur Conan Doyle'
		if corpus == 'austen':
			title = 'Pride and Prejudice'
		if corpus == 'bible':
			title = 'The Bible'
		if corpus == 'paper1':
			title =  str('PMID: '+ str(eligible_papers[0][1]))
		if corpus == 'paper2':
			title = str('PMID: '+ str(eligible_papers[1][1]))
		if corpus == 'paper3':
			title = str('PMID: ' + str(eligible_papers[2][1]))
		if corpus == 'paper4':
			title = str('PMID: ' + str(eligible_papers[3][1]))
		if corpus == 'paper5':
			title = str('PMID: ' + str(eligible_papers[4][1]))
		x, y, names, color = vis_scifi(corpus, query, eligible_papers)
		return render_template('results_scifi.html', x=x, y=y, title=title, color=color, query=query, names=names, eligible_papers=eligible_papers)
	else:
		logging.info("scifi analysis")
		corpus = 'darwin'
		title = 'On The Origin of Species'
		x, y, names, color = vis_scifi(corpus, query, eligible_papers)
		logging.info("done with x and y")
		if len(eligible_papers) < len(pmid_list):
			flash('Some input paper(s) are not avaliable for TextCompare')
		return render_template('results_scifi.html', x=x, y=y, title=title, color=color, query=query, names=names, eligible_papers=eligible_papers)

#################### OTHER ####################################################
@app.route('/testingstuff/')
def testingstuff():
	return render_template('test.html')


#Handles 404 errors
@app.errorhandler(404)
def page_not_found(e):
	return("you shouldnt be here!!")


#Configuration settings
if __name__ == '__main__':
	run_simple('0.0.0.0', 5000, app, use_reloader=True)
	#app.run(host='0.0.0.0') #dont want app.run() for uwsgi


#TODO: Add unit tests and such for Git & Travis UI
########### GRAVEYARD ##########################################################


