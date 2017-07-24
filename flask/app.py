from flask import Flask, render_template, request, flash, url_for, redirect, session, g, Blueprint
from flask_wtf import Form
from wtforms import TextField, SelectField
import gc, time, datetime, pickle, os.path, json
import sys, csv
from werkzeug.serving import run_simple
sys.path.append('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/')
from configapp import app, engine, connection, inputPapers
from content_management import * #mine
from citation_venn import make_venn #mine
from processors import *
from cache_lemma_nes import print_lemma_nes_samples, concat_lemma_nes_samples, exists_lemma, exists_nes


#Create Form for handling user-entered pmid
#Need to pass pmid in form to Entrez_IR.py
class pmidForm(Form):
	pmid = TextField('PubmedID')


#Main page
#Prints sample results from 2 coge publications
#User inputs a pubmed id and is then redirected to /results
@app.route("/cogecrawl/")
def cogecrawl():
	query = '18952863+18269575'
	coge_conn = connection()
	citations_with_links = db_unique_citations_retrieval(query, coge_conn) #unique
	unique_publications = db_unique_citations_number(query, coge_conn)
	coge_conn.close()
	return render_template("dashboard.html", citations_with_links=citations_with_links, unique_publications=unique_publications)


#Getting Flask-WTFs to work with sqlite3 here
#This function uses user entry to run the Entrez_IR.py 
#User entered pmid is entered into sqlite3 database
#TODO: make the results/ url unique 
@app.route('/results/', methods=["POST"])
def results():
	logging.info("In app route RESULTS")
	form = pmidForm()
	r_conn = connection() #results_connection to db
	try:
		if request.method == 'POST':
			entry = form.pmid.data #THIS IS THE USER INPUT FROM THE FORM #referencing 'class pmidForm'
			pmid_list = multiple_pmid_input(entry) #list for handling multiple pmids
			logging.info(pmid_list)


			# If the user inputs more than 5 PMIDs, return the home page and flash a warning
			# Need 5 PMIDs or less
			if len(pmid_list) > 5:
				flash('You have entered more than 5 PMIDs. Please reduce your query to 5 PMIDs or less to continue.')
				citations_with_links = db_unique_citations_retrieval(query, r_conn)  # unique
				unique_publications = db_unique_citations_number(query, r_conn)
				conn.close()
				return render_template("dashboard.html", citations_with_links=citations_with_links,
									   unique_publications=unique_publications)


			q = '+'
			query = str(q.join(pmid_list))
			logging.info("query: " + str(query))

			needed_to_annotate_check = []

			for user_input in pmid_list:
				logging.info(str(user_input))
				user_input = str(user_input)

				############################################
				#Check database for pmid #Does the entry exists in the db already?
				#conn = connection()
				s = inputPapers.select().\
					where(inputPapers.c.pmid == user_input)
				c = r_conn.execute(s)
				check1 = c.fetchone()
				c.close()


				#if the entry does NOT exist in the db already, will need to retrieve text, annotate it, and populate cache
				if check1 is None:
					update_check = "yes"
					flash('new pubmedid!')

					#Information Retireval of citing pmcids and info about them
					run_IR_not_db(user_input, r_conn)

					#Annotate
					logging.info("beginning multi-preprocessing")
					biodoc_data = do_multi_preprocessing(user_input, r_conn)
					needed_to_annotate_check_to_annotate.append("yes")
					logging.info("done with new document multi_preprocessing")
					logging.info("writing the BIODOC LEMMAS")

					#Populate cache (lemmas and nes)
					need_to_annotate = "yes" #of course we need to annotate, its a new pmid!
					print_lemma_nes_samples(user_input, biodoc_data, need_to_annotate)
					logging.info("* wrote lemme and nes samples to cache!!!")


					#After all citations have been processed, now we can do the analyses:
					if user_input == pmid_list[-1]: #if its the last pmid
						logging.info("last pmid in the query")

						#Lemma_samples and nes_samples for entire query here:
						logging.info("concatting lemma nes samples for query")
						#Populate cache
						need_to_update = "yes" #of course we do, something is new!
						concat_lemma_nes_samples(query, need_to_update)

						#Update "queries" table of db here!!
						logging.info("STARTING db_query_update_statistics")
						db_query_update_statistics(query, r_conn)
						logging.info("finishe db_query_update_statistics")


				#if the entry IS in the db, no need to retrieve text from Entrez, just grab from db
				#MAYBE need to annotate some new documents, maybe not
				#If new citations do need to be retireved, annotated, etc then DO NEED to re-populate cache
				if check1 is not None:
					update_check = "no" #no by default


					flash("alreay exists in database :) ")
					#Using user_input for Information Retireval - checks if any new papers have been added that we need to scrape
					need_to_annotate = run_IR_in_db(user_input, r_conn)


					if need_to_annotate == 'yes':
						needed_to_annotate_check.append('yes')
						logging.info("need to annotate new documents")
						#Annotate
						biodoc_data = do_multi_preprocessing(user_input, r_conn)
						logging.info("done with new document multi_preprocessing")
						#If need_to_annotate is "yes", will re-populate :)
						#Make cache
						print_lemma_nes_samples(user_input, biodoc_data, need_to_annotate)
						logging.info("repopulated lemmas and nes cache")

					#Make sure that the lemma and nes cache exists before moving on!!!
					if need_to_annotate == 'no':
						logging.info("dont need to annotate any new documents")
						if exists_lemma(user_input) and exists_nes(user_input):
							logging.info("lemmas and nes cache exist so pass :)")
							needed_to_annotate_check.append('no')
						else:
							logging.info("lemmas and nes cache didn't exist so gotta make them!!!")
							biodoc_data = do_multi_preprocessing(user_input, r_conn)
							logging.info("done with new document multi_preprocessing")
							nonexistant = "yes"
							print_lemma_nes_samples(user_input, biodoc_data, nonexistant)
							needed_to_annotate_check.append('yes')


					## Now that we have all the data, do the topic model
					## Only want to save final topic model (not running topic model)
					if user_input == pmid_list[-1]:
						logging.info("last pmid in the query")

						#If ANY user_inputs in the query needed to update, we must update the query's comprehensive cache.
						if 'yes' in needed_to_annotate_check:
							need_to_update = 'yes'
							#will over-ride existing file :)
							concat_lemma_nes_samples(query, need_to_update)
							update_check = "yes"

						if 'yes' not in needed_to_annotate_check:
							need_to_update = 'no'
							concat_lemma_nes_samples(query, need_to_update) #this should check for file
							update_check = "no"

						# Update "queries" table of db here!!
						logging.info("STARTING db_query_update_statistics")
						db_query_update_statistics(query, r_conn)
						logging.info("finishe db_query_update_statistics")


				#Housekeeping
				gc.collect() #garbage collector for cleaning up unneeded stuff
				session['entered_id'] = True
				session['engaged'] = 'engaged'


		citations_with_links = db_unique_citations_retrieval(query, r_conn) #unique
		# :( this s failing!!!!
		#unique_publications = db_unique_citations_number(query) <-- not sure why this is here...
		return render_template('results.html', form=form, citations_with_links=citations_with_links,
	   			query=query, update_check=update_check)

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

#NB: cogeembeddings does not connect to db
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

		filepath = embedding_lookup(query, k_clusters, window)
		logging.info(filepath)
		return render_template('coge_embeddings.html', filepath=filepath)
	else:
		filepath = 'coge_embed.json'
		return render_template('coge_embeddings.html', filepath=filepath)


#NB: cogelsa does not connect to db
#TODO: new default
@app.route('/cogelsa/', methods=["GET","POST"]) #default coge lsa for iframe
def cogelsa():
	form = visOptions()

	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		query = '18952863+18269575'
		k = int(k_clusters)

		try:
			jsonLSA= load_lsa(query, k)
		except Exception as e:
			# load data for analysis
			lemma_samples = load_lemma_cache(query)
			lsa_lemmas = [l[1] for l in lemma_samples]

			num_pubs = len(lsa_lemmas)
			if num_pubs < k:
				logging.info("k value is larger than number of publications")

			# run analysis & save it
			jsonLSA = run_lsa1(lsa_lemmas, k)
			print_lsa(query, jsonLSA, k)
		return render_template('coge_lsa.html', form=form, jsonLSA=jsonLSA)
	else:
		filename = "coge_lsa.json"
		filepath = os.path.join((app.config['PATH_TO_STATIC']), filename)
		with open(filepath) as load_data:
			jsonLSA = json.load(load_data)
		jsonLSA = re.sub('\'', '\"', str(jsonLSA)) #when I update the default, i can delete this trash!!
		return render_template('coge_lsa.html', form=form, jsonLSA=jsonLSA)


#NB: cogelda does not connect to db
#TODO: new default coge_lda
@app.route('/cogelda/', methods=["GET","POST"]) #default coge lda for iframe
def cogelda():
	form = visOptions()
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		num_words = form.w_words.data
		logging.info("the w value is "+str(num_words))
		query = '18952863+18269575'
		k = int(k_clusters)
		w = int(num_words)

		try:
			jsonLDA = load_lda(query, k, w)

		except Exception as e:
			lemma_samples = load_lemma_cache(query)
			lda_lemmas = [l[1] for l in lemma_samples]
			jsonLDA = run_lda1(lda_lemmas, k, w)
			print_lda(query, jsonLDA, k, w)

		return render_template('coge_lda.html', form=form, jsonLDA=jsonLDA)
	else:
		filename = "coge_lda1.json"
		filepath = os.path.join((app.config['PATH_TO_STATIC']), filename)
		with open(filepath) as load_data:
			jsonLDA = json.load(load_data)
		jsonLDA = re.sub('\'', '\"', str(jsonLDA)) #when I update the default, i can delete this trash!!
		return render_template('coge_lda.html', form=form, jsonLDA=jsonLDA)


#NBL cogejournals connects to db
@app.route('/cogejournals/') #default coge journals for iframe
def cogejournals():
	j_conn = connection()
	filename = "189/528/journals_18952863+18269575.json"
	filepath = os.path.join((app.config['PATH_TO_JOURNALS']), filename)
	with open(filepath) as load_data:
		journals = json.load(load_data)
	query = '18952863+18269575'
	range_years, unique_pubs, unique_journals = getJournalsVis(query, j_conn)
	years_list = range_years.split('+')
	s_year = years_list[0]
	e_year = years_list[1]
	j_conn.close()
	return render_template('coge_journals.html', journals=journals, unique_pubs=unique_pubs,
						   unique_journals=unique_journals, s_year=s_year, e_year=e_year)

#NB: wordcloud does NOT need db connection
@app.route('/cogewordcloud/', methods=["GET","POST"]) #default coge NES Word Cloud for iframe
def cogewordcloud():
	form = nesOptions()
	if request.method == 'POST':

		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		filename = "189/528/nes_18952863+18269575.pickle"
		nes_file = os.path.join((app.config['PATH_TO_CACHE']), filename)
		with open(nes_file, "rb") as f:
			nes_samples = pickle.load(f)

		nes_list = [n[1] for n in nes_samples]
		#id_list = [n[0] for n in nes_samples]

		wordcloud_data = vis_wordcloud(nes_list, nes_categories, w_number)
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong>Displaying results for N= '+str(w_number)+' from categories: '+str(nes_categories)+'</div>'
		return render_template('coge_wordcloud.html',  wordcloud_data=wordcloud_data, popup=popup)
	else:
		#Default data
		filename = "coge_wcloud1.json"
		filepath = os.path.join((app.config['PATH_TO_STATIC']), filename)
		with open(filepath) as load_data:
			wordcloud_data = json.load(load_data)
		wordcloud_data = re.sub('\'', '\"', str(wordcloud_data)) #for some reason JS wants the \ by the quotes...
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong> Default: N=10, from all categories.</div>'
		return render_template('coge_wordcloud.html', wordcloud_data=wordcloud_data, popup=popup)

#NB: uses db
@app.route('/cogeheatmap/', methods=["GET","POST"]) #default coge NES heatmap for iframe
def cogeheatmap():
	form = nesOptions()
	query = '18952863+18269575'
	if request.method == 'POST':
		hm_conn = connection()
		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		filename1 = "189/528/nes_18952863+18269575.pickle"
		nes_file = os.path.join((app.config['PATH_TO_CACHE']), filename1)
		with open(nes_file, "rb") as f:
			nes_samples = pickle.load(f)

		filename2 = "189/528/lemma_samples_18952863+18269575.pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename2)
		with open(lemma_file, "rb") as f:
			lemma_samples = pickle.load(f)

		x_docs, y_words, z_counts, titles = vis_heatmap(lemma_samples, nes_samples, nes_categories, w_number, hm_conn)
		len_x_docs = list(range(len(x_docs)))
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong>Displaying results for N= '+str(w_number)+' from categories: '+str(nes_categories)+'</div>'
		hm_conn.close()
		return render_template('coge_heatmap2.html', z_counts=z_counts, x_docs=x_docs, y_words=y_words, titles=titles, len_x_docs=len_x_docs, popup=popup)
	else:
		#display the default data :D
		return render_template('coge_heatmap1.html')

#NB: uses db
@app.route('/cogeclustermap/', methods=["GET","POST"]) #default coge clustermap
def cogeclustermap():
	query = '18952863+18269575'
	form = nesOptions()
	if request.method == 'POST':
		cm_conn = connection()
		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is " + str(w_number))

		filename1 = "189/528/nes_18952863+18269575.pickle"
		nes_file = os.path.join((app.config['PATH_TO_CACHE']), filename1)
		with open(nes_file, "rb") as f:
			nes_samples = pickle.load(f)

		filename2 = "189/528/lemma_samples_18952863+18269575.pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename2)
		with open(lemma_file, "rb") as f:
			lemma_samples = pickle.load(f)

		saveName = vis_clustermap(lemma_samples, nes_samples, nes_categories, w_number, query, cm_conn)
		#Fix.... why bother returning full file path if html only looks in "static" dir?
		image = "clustermaps/cm_"+str(query)+".png"
		cm_conn.close()
		return render_template('coge_clustermap.html', image=image)
	else:
		#show the default data :D
		image = "clustermaps/cm_18952863+18269575.png"
		return render_template('coge_clustermap.html', image=image)

# NB: uses db
@app.route('/cogekmeans/', methods=["GET","POST"]) #default coge k-means clustering for iframe
def cogekmeans():
	form = visOptions()
	if request.method == 'POST':
		k_conn = connection()
		filename = "189/528/lemma_samples_18952863+18269575.pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename)
		with open(lemma_file, "rb") as f:
			lemma_samples = pickle.load(f)

		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates, titles0, titles1, titles2, titles3, titles4 = vis_kmeans(lemma_samples, k_clusters, k_conn)
		k_conn.close()
		return render_template('coge_kmeans2.html', x0_coordinates=x0_coordinates, y0_coordinates=y0_coordinates, z0_coordinates=z0_coordinates,
							   x1_coordinates=x1_coordinates, y1_coordinates=y1_coordinates, z1_coordinates=z1_coordinates,
							   x2_coordinates=x2_coordinates, y2_coordinates=y2_coordinates, z2_coordinates=z2_coordinates,
							   x3_coordinates=x3_coordinates, y3_coordinates=y3_coordinates, z3_coordinates=z3_coordinates,
							   x4_coordinates=x4_coordinates, y4_coordinates=y4_coordinates, z4_coordinates=z4_coordinates,
							   titles0=titles0, titles1=titles1, titles2=titles2, titles3=titles3, titles4=titles4)
	else:
		#show the default data!!
		return render_template('coge_kmeans.html')


#NB: uses db
@app.route('/coge_stats/') #default coge statistics for iframe
def coge_stats():
	s_conn = connection()
	query = "18952863+18269575"
	input_click_citations = statsSelfInfo(query, s_conn)
	statistics = get_statistics(query, s_conn) #actually a lot of these are "None" right now. Will need to populate.
	sum_total = statistics[0]
	unique = statistics[1]
	sum_abstracts = statistics[2]
	sum_whole = statistics[3]
	sum_sents = statistics[4]
	sum_tokens = statistics[5]
	s_conn.close()
	return render_template('coge_stats.html', input_click_citations=input_click_citations,
						   sum_total=sum_total, unique=unique, sum_abstracts=sum_abstracts, sum_whole=sum_whole,
						   sum_sents=sum_sents, sum_tokens=sum_tokens)


#NB: uses db
@app.route('/coge_scifi/', methods=["GET","POST"]) #default coge scifi for iframe
def coge_scifi():
	form = corpusOptions()
	path_to_eligible_paper = os.path.join((app.config['PATH_TO_CACHE']), '259/367/2593677.txt')
	eligible_papers = [('paper1', '18952863', path_to_eligible_paper, '2008, Lyons')]
	if request.method == 'POST':
		csf_conn = connection()
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
			title = 'The Brain & The Voice in Speech & Song by Mott Frederick Walker'
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
			title =  'Lyons et al., 2008'
		x, y, names, color = vis_scifi(corpus, query, eligible_papers, csf_conn)
		csf_conn.close()
		return render_template('coge_scifi2.html', x=x, y=y, title=title, color=color, names=names, eligible_papers=eligible_papers)
	else:
		flash('Some input paper(s) are not avaliable for TextCompare')
		#Default data!
		return render_template('coge_scifi.html', eligible_papers=eligible_papers)


############### Results visualizations #########################################


#NB: needs db connection
@app.route('/resjournals/<query>/<update_check>', methods=["GET", "POST"]) #user journals for iframe
def resjournals(query, update_check):
	#need to get last user_input
	logging.info("in routine res-journals")
	jr_conn = connection() #journal results connection

	needed_to_annotate_check = [update_check]
	range_years, unique_publications, unique_journals = print_journalvis(query, needed_to_annotate_check, jr_conn)

	logging.info("YEARS RANGE: " +str(range_years))
	#Need years for range
	### AHHHH if range is not in db yet, it will get it without the "+" as a tuple!
	try:
		years_list = range_years.split('+')
		s_year = years_list[0]
		e_year = years_list[1]
	except Exception as e:
		s_year = years_list[0]
		e_year = years_list[1]

	#only want to load the json for the LAST id in the query (so includes all)
	pmid_list = query.split('+') #list of string pmids
	pmid = pmid_list[0]
	prefix = pmid[0:3]
	suffix = pmid[3:6]

	filename = str(prefix)+'/'+str(suffix)+'/'+"journals_"+str(query)+".json"
	logging.info("last entry's JOURNAL is named: " + str(filename))
	savePath = (app.config['PATH_TO_JOURNALS'])
	completeName = os.path.join(savePath, filename)
	logging.info("complete file: " + str(completeName))
	with open(completeName) as load_data:
		journals = json.load(load_data)
	jr_conn.close()
	return render_template('results_journals.html', journals=journals, s_year=s_year, e_year=e_year,
						   unique_journals=unique_journals, unique_publications=unique_publications)



#NB: does NOT need a db connection
@app.route('/resembed/<query>/<update_check>', methods=["GET", "POST"]) #user embeddings for iframe
def resembeddings(query, update_check):
	form = visOptions()
	if request.method == 'POST':
		window = int(form.w_words.data)
		logging.info(window)
		k_clusters = int(form.k_val.data)  # 2,3,4,or 5
		logging.info(k_clusters)

		pmid_list = query.split('+')  # list of string pmids
		pmid = pmid_list[0]  # get the first
		prefix = pmid[0:3]
		suffix = pmid[3:6]

		if update_check == 'yes':
			run_embeddings(query, k_clusters, window)  # 50 words in 6 clusters
			filename = str(prefix) + '/' + str(suffix) + '/' + 'fgraph_' + str(query) + '_' + str(
				k_clusters) + '_' + str(window) + '.json'
			filepath = os.path.join((app.config['PATH_TO_FGRAPHS']), filename)

		if update_check == 'no':
			filepath = embedding_lookup(query, k_clusters, window)

		return render_template('results_embeddings.html', query=query, update_check=update_check, filepath=filepath)
	else:
		k_clusters = 10
		window = 100
		logging.info("update_check: " + str(update_check))

		pmid_list = query.split('+')  # list of string pmids
		pmid = pmid_list[0]  # get the first
		prefix = pmid[0:3]
		suffix = pmid[3:6]

		if update_check == 'yes':
			run_embeddings(query, k_clusters, window)  # 50 words in 6 clusters
			filename = str(prefix) + '/' + str(suffix) + '/' + 'fgraph_' + str(query) + '_' + str(
				k_clusters) + '_' + str(window) + '.json'
			filepath = os.path.join('fgraphs', filename)

		if update_check == 'no':
			filepath = embedding_lookup(query, k_clusters, window)

		return render_template('results_embeddings.html', query=query, update_check=update_check, filepath=filepath)


#NB: does NOT need a db connection
@app.route('/reslsa/<query>/<update_check>', methods=["GET", "POST"]) #user lsa for iframe
def reslsa(query, update_check):
	form = visOptions()
	if request.method == 'POST':
		k_clusters = form.k_val.data  # 2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		k = int(k_clusters)

		#update check? yes --> force update of file and save it
		#update check? no --> look for file, if no file, make it
		if update_check == 'no':
			#look for file, if no file, make it
			jsonDict = load_lsa(query, k)

		if update_check == 'yes':
			logging.info("LSA: RERUN")

			#load data for analysis
			lemma_samples = load_lemma_cache(query)
			lsa_lemmas = [l[1] for l in lemma_samples]

			num_pubs = len(lsa_lemmas)
			if num_pubs < k:
				logging.info("k value is larger than number of publications")

			#run analysis
			jsonDict = run_lsa1(lsa_lemmas, k)
			#save it
			print_lsa(query, jsonDict, k)

			logging.info("did it all!")
		return render_template('results_lsa.html', query=query, jsonDict=jsonDict, update_check=update_check)
	else:
		logging.info("LSA: DEFAULT (results)")
		logging.info("LSA ID: " +str(query))

		if update_check == 'yes':
			#force update
			##load the lemmas
			k=7
			lemma_samples = load_lemma_cache(query)

			lsa_lemmas = [l[1] for l in lemma_samples]

			jsonDict = run_lsa1(lsa_lemmas, k)
			print_lsa(query, jsonDict, k)


		if update_check == 'no':
			k=7
			jsonDict = load_lsa(query, k)

		return render_template('results_lsa.html', query=query, jsonDict=jsonDict, update_check=update_check)


#NB: does NOT need a db connection
@app.route('/reslda/<query>/<update_check>', methods=["GET", "POST"]) #user lda for iframe
def reslda(query, update_check):
	form = visOptions()
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		num_words = form.w_words.data
		logging.info("the w value is "+str(num_words))
		k = int(k_clusters)
		w = int(num_words)

		if update_check == 'yes':
			lemma_samples = load_lemma_cache(query)
			lda_lemmas = [l[1] for l in lemma_samples]
			jsonLDA = run_lda1(lda_lemmas, k, w)
			print_lda(query, jsonLDA, k, w)


		if update_check == 'no':
			jsonLDA = load_lda(query, k, w)


		return render_template('results_lda.html', form=form, jsonLDA=jsonLDA, query=query, update_check=update_check)
	else:
		#need to get last user_input
		#use id to do stuff
		logging.info("in routine reslDa")
		logging.info("RES-LDA ID: " +str(query))
		k = 7
		w = 7

		if update_check == 'yes':
			lemma_samples = load_lemma_cache(query)
			lda_lemmas = [l[1] for l in lemma_samples]
			jsonLDA = run_lda1(lda_lemmas, k, w)
			print_lda(query, jsonLDA, k, w)


		if update_check == 'no':
			jsonLDA = load_lda(query, k, w)

		return render_template('results_lda.html', form=form, jsonLDA=jsonLDA, query=query, update_check=update_check)


#NB: does NOT need a db connection
@app.route('/reswordcloud/<query>', methods=["GET", "POST"]) #user wordcloud for iframe
def reswordcloud(query):
	form = nesOptions()
	if request.method == 'POST':

		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		pmid_list = query.split('+')
		pmid = pmid_list[0]
		prefix = pmid[0:3]
		suffix = pmid[3:6]

		filename1 = str(prefix) + '/' + str(suffix) + '/' + "nes_" + str(query) + ".pickle"
		nes_file = os.path.join((app.config['PATH_TO_CACHE']), filename1)
		with open(nes_file, "rb") as f:
			nes_samples = pickle.load(f)

		nes_list = [n[1] for n in nes_samples]
		wordcloud_data = vis_wordcloud(nes_list, nes_categories, w_number)
		popup = ' '
		return render_template('results_wordcloud.html', query=query,  wordcloud_data=wordcloud_data, popup=popup)
	else:
		nes_categories= ['BioProcess', 'CellLine', 'Cellular_component', 'Family', 'Gene_or_gene_product', 'Organ', 'Simple_chemical', 'Site', 'Species', 'TissueType']
		logging.info(nes_categories)
		w_number = 10
		logging.info("the w value is "+str(w_number))

		pmid_list = query.split('+')
		pmid = pmid_list[0]
		prefix = pmid[0:3]
		suffix = pmid[3:6]

		filename1 = str(prefix) + '/' + str(suffix) + '/' + "nes_" + str(query) + ".pickle"
		nes_file = os.path.join((app.config['PATH_TO_CACHE']), filename1)
		with open(nes_file, "rb") as f:
			nes_samples = pickle.load(f)

		nes_list = [n[1] for n in nes_samples]
		wordcloud_data = vis_wordcloud(nes_list, nes_categories, w_number)
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong> Default: N=10, from all categories.</div>'

		return render_template('results_wordcloud.html', query=query, wordcloud_data=wordcloud_data, popup=popup)


#NB: needs a db connection
#TODO: I have no idea where/when to close the connection
@app.route('/res_heatmap/<query>', methods=["GET", "POST"]) #user heatmap for iframe
def res_heatmap(query):
	form = nesOptions()
	hmr_conn = connection() #heatmap results database connection
	if request.method == 'POST':
		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		pmid_list = query.split('+')
		pmid = pmid_list[0]
		prefix = pmid[0:3]
		suffix = pmid[3:6]

		filename1 = str(prefix) + '/' + str(suffix) + '/' + "nes_" + str(query) + ".pickle"
		nes_file = os.path.join((app.config['PATH_TO_CACHE']), filename1)
		with open(nes_file, "rb") as f:
			nes_samples = pickle.load(f)

		filename2 = str(prefix) + '/' + str(suffix) + '/' + "lemma_samples_" + str(query) + ".pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename2)
		with open(lemma_file, "rb") as l:
			lemma_samples = pickle.load(l)

		x_docs, y_words, z_counts, titles = vis_heatmap(lemma_samples, nes_samples, nes_categories, w_number, hmr_conn)
		popup = ' '
		hmr_conn.close()
		return render_template('results_heatmap.html', query=query, z_counts=z_counts, x_docs=x_docs, y_words=y_words, popup=popup, titles=titles)
	else:
		nes_categories= ['BioProcess', 'CellLine', 'Cellular_component', 'Family', 'Gene_or_gene_product', 'Organ', 'Simple_chemical', 'Site', 'Species', 'TissueType']
		w_number = 10
		logging.info("the w value is "+str(w_number))

		pmid_list = query.split('+')
		pmid = pmid_list[0]
		prefix = pmid[0:3]
		suffix = pmid[3:6]

		filename1 = str(prefix) + '/' + str(suffix) + '/' + "nes_" + str(query) + ".pickle"
		nes_file = os.path.join((app.config['PATH_TO_CACHE']), filename1)
		with open(nes_file, "rb") as f:
			nes_samples = pickle.load(f)

		filename2 = str(prefix) + '/' + str(suffix) + '/' + "lemma_samples_" + str(query) + ".pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename2)
		with open(lemma_file, "rb") as l:
			lemma_samples = pickle.load(l)

		x_docs, y_words, z_counts, titles = vis_heatmap(lemma_samples, nes_samples, nes_categories, w_number, hmr_conn)
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong> Default: N=10, from all categories.</div>'
		hmr_conn.close()
		return render_template('results_heatmap.html', query=query, z_counts=z_counts, x_docs=x_docs, y_words=y_words, popup=popup, titles=titles)


#NB: needs a db connection
@app.route('/res_clustermap/<query>', methods=["GET", "POST"]) #user heatmap for iframe
def res_clustermap(query):
	form = nesOptions()
	if request.method == 'POST':
		cmr_conn = connection() #clustermap results database connection
		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		pmid_list = query.split('+')
		pmid = pmid_list[0]
		prefix = pmid[0:3]
		suffix = pmid[3:6]

		filename1 = str(prefix) + '/' + str(suffix) + '/' + "nes_" + str(query) + ".pickle"
		nes_file = os.path.join((app.config['PATH_TO_CACHE']), filename1)
		with open(nes_file, "rb") as f:
			nes_samples = pickle.load(f)

		filename2 = str(prefix) + '/' + str(suffix) + '/' + "lemma_samples_" + str(query) + ".pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename2)
		with open(lemma_file, "rb") as l:
			lemma_samples = pickle.load(l)

		saveName = vis_clustermap(lemma_samples, nes_samples, nes_categories, w_number, query, cmr_conn)
		image = '/clustermaps/' + saveName
		cmr_conn.close()
		return render_template('results_clustermap.html',image=image, query=query)
	else:
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong> Choose N and categories to run clustermap.</div>'
		#No default visualization
		return render_template('results_clustermapH.html', query=query, popup=popup)



#NB: needs a db connection
@app.route('/res_kmeans/<query>', methods=["GET", "POST"]) #user k-means for iframe
def res_kmeans(query):
	form = visOptions()
	if request.method == 'POST':
		kmr_conn = connection() #k-means results database connection
		pmid_list = query.split('+') #list of string pmids
		pmid = pmid_list[0]
		prefix = pmid[0:3]
		suffix = pmid[3:6]

		filename2 = str(prefix) + '/' + str(suffix) + '/' + "lemma_samples_" + str(query) + ".pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename2)
		with open(lemma_file, "rb") as l:
			lemma_samples = pickle.load(l)

		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates, titles0, titles1, titles2, titles3, titles4 = vis_kmeans(lemma_samples, k_clusters, kmr_conn)
		kmr_conn.close()
		return render_template('res_kmeans1.html', query=query,
		   x0_coordinates=x0_coordinates, y0_coordinates=y0_coordinates, z0_coordinates=z0_coordinates,
		   x1_coordinates=x1_coordinates, y1_coordinates=y1_coordinates, z1_coordinates=z1_coordinates,
		   x2_coordinates=x2_coordinates, y2_coordinates=y2_coordinates, z2_coordinates=z2_coordinates,
		   x3_coordinates=x3_coordinates, y3_coordinates=y3_coordinates, z3_coordinates=z3_coordinates,
		   x4_coordinates=x4_coordinates, y4_coordinates=y4_coordinates, z4_coordinates=z4_coordinates,
			titles0=titles0, titles1=titles1, titles2=titles2, titles3=titles3, titles4=titles4)
	else:
		#Not running any default vis because its slow
		return render_template('res_kmeans1.html', query=query)



@app.route('/res_stats/<query>', methods=["GET"]) #user statistics for iframe
def res_stats(query):
	sr_conn = connection() #statistics results db connection
	pmid_list = query.split('+') #list of string pmids
	venn_data = make_venn(pmid_list, sr_conn)

	input_click_citations = statsSelfInfo(query, sr_conn)
	statistics = get_statistics(query, sr_conn)
	sum_total = statistics[0]
	unique = statistics[1]
	sum_abstracts = statistics[2]
	sum_whole = statistics[3]
	sum_sents = statistics[4]
	sum_tokens = statistics[5]

	#get x, y coordinates for pubs x year bar chart.
	#max 5 papers
	x0, x1, x2, x3, x4, y0, y1, y2, y3, y4, n0, n1, n2, n3, n4 = stats_barchart(query, sr_conn)
	sr_conn.close()
	return render_template('results_stats.html', input_click_citations=input_click_citations,
						   venn_data=venn_data, sum_total=sum_total,
						   unique=unique, sum_abstracts=sum_abstracts, sum_whole=sum_whole,
						   sum_sents=sum_sents, sum_tokens=sum_tokens,
						   x0=x0, x1=x1, x2=x2, x3=x3, x4=x4,
						   y0=y0, y1=y1, y2=y2, y3=y3, y4=y4,
						   n0=n0, n1=n1, n2=n2, n3=n3, n4=n4)


#NB: needs db connection
#TODO: Idk when to close connection
@app.route('/results_scifi/<query>', methods=["GET","POST"]) #default coge scifi for iframe
def results_scifi(query):
	form = corpusOptions()
	pmid_list = query.split('+')  # list of string pmids
	#decide eligible papers:
	sfr_conn = connection() #scifi results db connection
	eligible_papers = inputEligible(query, sfr_conn)

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
		x, y, names, color = vis_scifi(corpus, query, eligible_papers, sfr_conn)
		sfr_conn.close()
		return render_template('results_scifi.html', x=x, y=y, title=title, color=color, query=query, names=names, eligible_papers=eligible_papers)
	else:
		logging.info("scifi analysis")
		corpus = 'darwin'
		title = 'On The Origin of Species'
		x, y, names, color = vis_scifi(corpus, query, eligible_papers, sfr_conn)
		logging.info("done with x and y")
		if len(eligible_papers) < len(pmid_list):
			flash('Some input paper(s) are not avaliable for TextCompare')
		sfr_conn.close()
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
	#app.run(host='0.0.0.0') if you are not running the app with uwsgi!


#TODO: Add unit tests and such for Git & Travis UI
########### GRAVEYARD ##########################################################


