from flask import Flask, render_template, request, flash, url_for, redirect, session, g, Blueprint
from flask_wtf import Form
from wtforms import TextField, SelectField
import gc, time, datetime, pickle, os.path, json
import sys, csv
from werkzeug.serving import run_simple
sys.path.append('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/')
from database_management import engine, connection, inputPapers #mine
from content_management import * #mine
from citation_venn import make_venn #mine
from processors import *
from cache_lemma_nes import print_lemma_nes_samples, concat_lemma_nes_samples


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
	query = '18952863+18269575'
	citations_with_links = db_unique_citations_retrieval(query) #unique
	unique_publications = db_unique_citations_number(query)
	return render_template("dashboard.html", citations_with_links=citations_with_links, unique_publications=unique_publications)


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


				#if the entry does NOT exist in the db already, will need to retrieve text, annotate it, and populate cache
				if check1 is None:
					update_check = "yes"
					flash('new pubmedid!')
					#Using user_input for Information Retireval of citing pmcids and info about them
					run_IR_not_db(user_input)
					logging.info("beginning multi-preprocessing")
					#Annotate
					biodoc_data = do_multi_preprocessing(user_input)
					logging.info("done with new document multi_preprocessing")
					logging.info("writing the BIODOC LEMMAS")
					#Populate cache (lemmas and nes)
					biodoc_to_db(biodoc_data) #writes data_samples (lemmas) and NER to db
					logging.info("* done writing the biodoc lemmas")
					logging.info("* writing to cache: lemma_samples and nes_samples ")
					print_lemma_nes_samples(user_input, biodoc_data)
					logging.info("* wrote lemme and nes samples to cache!!!")


					#After all citations have been processed, now we can do the analyses:
					if user_input == pmid_list[-1]: #if its the last pmid
						logging.info("last pmid in the query")

						#Lemma_samples and nes_samples for entire query here:
						logging.info("concatting lemma nes samples for query")
						#Populate cache
						concat_lemma_nes_samples(query)

						#Update "queries" table of db here!!
						db_query_update_statistics(query)

						# logging.info(user_input+" is the last one (LDA)")
						# jsonLDA = run_lda1(data_samples, 3, 5)
						# print_lda(query, user_input, jsonLDA) #print lda topic model to json



				#if the entry IS in the db, no need to retrieve text from Entrez, just grab from db
				#MAYBE need to annotate some new documents, maybe not
				#If new citations do need to be retireved, annotated, etc then DO NEED to re-populate cache
				if check1 is not None:
					update_check = "no" #no by default

					needed_to_annotate_check = []

					flash("alreay exists in database :) ")
					#Using user_input for Information Retireval - checks if any new papers have been added that we need to scrape
					need_to_annotate = run_IR_in_db(user_input)

					if need_to_annotate == 'yes':
						needed_to_annotate_check.append('yes')
						logging.info("need to annotate new documents")
						biodoc_data = do_multi_preprocessing(user_input)
						logging.info("done with new document multi_preprocessing")
						#If need_to_annotate is "yes", will re-populate :)
						print_lemma_nes_samples(user_input, biodoc_data, need_to_annotate)
						logging.info("repopulated lemmas and nes cache")

					#Don't need to re-populat cache
					#TODO: maybe still double check if the lemma & nes files exist?? :s
					if need_to_annotate == 'no':
						logging.info("dont need to annotate any new documents")
						needed_to_annotate_check.append('no')
						pass


					## Now that we have all the data, do the topic model
					## Only want to save final topic model (not running topic model)
					if user_input == pmid_list[-1]:
						logging.info("last pmid in the query")

						#If ANY user_inputs in the query needed to update, we must update the query's comprehensive cache.
						if 'yes' in needed_to_annotate_check:
							need_to_update = 'yes'
							#will over-ride existing file :)
							concat_lemma_nes_samples(query, need_to_update)
							db_query_update_statistics(query)
							update_check = "yes"

						if 'yes' not in needed_to_annotate_check:
							#update_check is "no" by default
							pass

                        #
						# logging.info(user_input + " is the last one (LDA)")
						# jsonLDA = run_lda1(data_samples, 3, 5)
						# print_lda(query, user_input, jsonLDA)  # print


				#Housekeeping
				gc.collect() #garbage collector for cleaning up unneeded stuff
				session['entered_id'] = True
				session['engaged'] = 'engaged'


		citations_with_links = db_unique_citations_retrieval(query) #unique
		unique_publications = db_unique_citations_number(query)
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

#TODO: just load cache
@app.route('/cogembeddings/', methods=["GET","POST"]) #default coge embeddings topic for iframe
def cogeembeddings():
	form = visOptions()
	if request.method =='POST':
		logging.info("posted something to cogembeddings")
		query = '18952863+18269575'
		# pmid_list = query.split('+')  # list of string pmids
		# pmid = pmid_list[0]  # get the first
		# prefix = pmid[0:3]
		# suffix = pmid[3:6]

		window = int(form.w_words.data)
		logging.info(window)
		k_clusters = int(form.k_val.data)  # 2,3,4,or 5
		logging.info(k_clusters)

		# run_embeddings(query, k_clusters, window)  # 50 words in 6 clusters
		# filename = str(prefix) + '/' + str(suffix) + '/' + 'fgraph_' + str(query) + '_' + str(k_clusters) + '_' + str(
		# 	window) + '.json'
        #
		# filepath = os.path.join('fgraphs', filename) #should this be os.path.abspath()?
		filepath = embedding_lookup(query, k_clusters, window)

		return render_template('coge_embeddings.html', filepath=filepath)
	else:
		filepath = 'coge_embed.json'
		return render_template('coge_embeddings.html', filepath=filepath)


#TODO: new default
#TODO: update for chaching sub options
@app.route('/cogelsa/', methods=["GET","POST"]) #default coge lsa for iframe
def cogelsa():
	form = visOptions()
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		query = '18952863+18269575'

		filename = '189/528/lemma_samples_18952863+18269575.pickle'
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename)
		with open(lemma_file, "rb") as f:
			lemma_samples = pickle.load(f)

		logging.info("GETTING THE DATA ....")
		lemmas_for_lsa = [l[1] for l in lemma_samples] #ignore the pmcid's in l[0], ignore tags in l[2]

		logging.info("rerunning the analysis")
		k = int(k_clusters)

		jsonLSA = run_lsa1(lemmas_for_lsa, k)
		logging.info("did it all!")
		return render_template('coge_lsa.html', form=form, jsonLSA=jsonLSA)
	else:
		filename = "coge_lsa.json"
		filepath = os.path.join((app.config['PATH_TO_STATIC']), filename)
		with open(filepath) as load_data:
			jsonLSA = json.load(load_data)
		jsonLSA = re.sub('\'', '\"', str(jsonLSA)) #when I update the default, i can delete this trash!!
		return render_template('coge_lsa.html', form=form, jsonLSA=jsonLSA)


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

		filename = "189/528/lemma_samples_18952863+18269575.pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename)
		with open(lemma_file, "rb") as f:
			lemma_samples = pickle.load(f)

		lemmas_for_lda = [l[1] for l in lemma_samples] #ignore the pmcid's in l[0], ignore tags in l[2]

		logging.info("rerunning the analysis")
		k = int(k_clusters)
		w = int(num_words)
		jsonLDA = run_lda1(lemmas_for_lda, k, w)
		return render_template('coge_lda.html', form=form, jsonLDA=jsonLDA)
	else:
		filename = "coge_lda1.json"
		filepath = os.path.join((app.config['PATH_TO_STATIC']), filename)
		with open(filepath) as load_data:
			jsonLDA = json.load(load_data)
		jsonLDA = re.sub('\'', '\"', str(jsonLDA)) #when I update the default, i can delete this trash!!
		return render_template('coge_lda.html', form=form, jsonLDA=jsonLDA)


@app.route('/cogejournals/') #default coge journals for iframe
def cogejournals():
	filename = "189/528/journals_18952863+18269575.json"
	filepath = os.path.join((app.config['PATH_TO_JOURNALS']), filename)
	with open(filepath) as load_data:
		journals = json.load(load_data)
	query = '18952863+18269575'
	range_years, unique_pubs, unique_journals = getJournalsVis(query)
	years_list = range_years.split('+')
	s_year = years_list[0]
	e_year = years_list[1]
	return render_template('coge_journals.html', journals=journals, unique_pubs=unique_pubs,
						   unique_journals=unique_journals, s_year=s_year, e_year=e_year)


@app.route('/cogewordcloud/', methods=["GET","POST"]) #default coge NES Word Cloud for iframe
def cogewordcloud():
	form = nesOptions()
	if request.method == 'POST':

		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		filename = "189/528/nes_18952863+18269575.pickle"
		#nes_file = '/home/hclent/data/pmcids/189/528/nes_18952863+18269575.pickle'
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


@app.route('/cogeheatmap/', methods=["GET","POST"]) #default coge NES heatmap for iframe
def cogeheatmap():
	form = nesOptions()
	query = '18952863+18269575'
	if request.method == 'POST':
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

		x_docs, y_words, z_counts, titles = vis_heatmap(lemma_samples, nes_samples, nes_categories, w_number)
		len_x_docs = list(range(len(x_docs)))
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong>Displaying results for N= '+str(w_number)+' from categories: '+str(nes_categories)+'</div>'
		return render_template('coge_heatmap2.html', z_counts=z_counts, x_docs=x_docs, y_words=y_words, titles=titles, len_x_docs=len_x_docs, popup=popup)
	else:
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

		filename1 = "189/528/nes_18952863+18269575.pickle"
		nes_file = os.path.join((app.config['PATH_TO_CACHE']), filename)
		with open(nes_file, "rb") as f:
			nes_samples = pickle.load(f)

		filename2 = "189/528/lemma_samples_18952863+18269575.pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename2)
		with open(lemma_file, "rb") as f:
			lemma_samples = pickle.load(f)

		saveName = vis_clustermap(lemma_samples, nes_samples, nes_categories, w_number, query)
		#Fix.... why bother returning full file path if html only looks in "static" dir?
		image = "clustermaps/cm_"+str(query)+".png"
		return render_template('coge_clustermap.html', image=image)
	else:
		image = "clustermaps/cm_18952863+18269575.png"
		return render_template('coge_clustermap.html', image=image)


@app.route('/cogekmeans/', methods=["GET","POST"]) #default coge k-means clustering for iframe
def cogekmeans():
	form = visOptions()
	if request.method == 'POST':

		filename = "189/528/lemma_samples_18952863+18269575.pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename)
		with open(lemma_file, "rb") as f:
			lemma_samples = pickle.load(f)

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


#TODO: code to auto fill in venn diagram?
@app.route('/coge_stats/') #default coge statistics for iframe
def coge_stats():
	query = "18952863+18269575"
	input_click_citations = statsSelfInfo(query)
	statistics = get_statistics(query) #actually a lot of these are "None" right now. Will need to populate.
	sum_total = statistics[0]
	unique = statistics[1]
	sum_abstracts = statistics[2]
	sum_whole = statistics[3]
	sum_sents = statistics[4]
	sum_tokens = statistics[5]
	return render_template('coge_stats.html', input_click_citations=input_click_citations,
						   sum_total=sum_total, unique=unique, sum_abstracts=sum_abstracts, sum_whole=sum_whole,
						   sum_sents=sum_sents, sum_tokens=sum_tokens)


@app.route('/coge_scifi/', methods=["GET","POST"]) #default coge scifi for iframe
def coge_scifi():
	form = corpusOptions()
	path_to_eligible_paper = os.path.join((app.config['PATH_TO_CACHE']), '259/367/2593677.txt')
	eligible_papers = [('paper1', '18952863', path_to_eligible_paper, '2008, Lyons')]
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
			title =  'Lyons et al., 2008'
		x, y, names, color = vis_scifi(corpus, query, eligible_papers)
		return render_template('coge_scifi2.html', x=x, y=y, title=title, color=color, names=names, eligible_papers=eligible_papers)
	else:
		flash('Some input paper(s) are not avaliable for TextCompare')
		return render_template('coge_scifi.html', eligible_papers=eligible_papers)


############### Results visualizations #########################################
@app.route('/resjournals/<query>/<update_check>', methods=["GET", "POST"]) #user journals for iframe
def resjournals(query, update_check):
	#need to get last user_input
	logging.info("in routine res-journals")

	needed_to_annotate_check = [update_check]
	range_years, unique_publications, unique_journals = print_journalvis(query, needed_to_annotate_check)

	logging.info("YEARS RANGE: " +str(range_years))
	#Need years for range
	years_list = range_years.split('+')
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
	return render_template('results_journals.html', journals=journals, s_year=s_year, e_year=e_year,
						   unique_journals=unique_journals, unique_publications=unique_publications)



#TODO: re-implement resembeddings for results!
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
		logging.info("update_check: ", update_check)

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


@app.route('/res_heatmap/<query>', methods=["GET", "POST"]) #user heatmap for iframe
def res_heatmap(query):
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

		filename2 = str(prefix) + '/' + str(suffix) + '/' + "lemma_samples_" + str(query) + ".pickle"
		lemma_file = os.path.join((app.config['PATH_TO_CACHE']), filename2)
		with open(lemma_file, "rb") as l:
			lemma_samples = pickle.load(l)

		x_docs, y_words, z_counts, titles = vis_heatmap(lemma_samples, nes_samples, nes_categories, w_number)
		popup = ' '
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

		x_docs, y_words, z_counts, titles = vis_heatmap(lemma_samples, nes_samples, nes_categories, w_number)
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

		saveName = vis_clustermap(lemma_samples, nes_samples, nes_categories, w_number, query)
		image = '/clustermaps/' + saveName

		return render_template('results_clustermap.html',image=image, query=query)
	else:
		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong> Choose N and categories to run clustermap.</div>'
		return render_template('results_clustermapH.html', query=query, popup=popup)


@app.route('/res_kmeans/<query>', methods=["GET", "POST"]) #user k-means for iframe
def res_kmeans(query):
	form = visOptions()
	if request.method == 'POST':

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
		x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates, titles0, titles1, titles2, titles3, titles4 = vis_kmeans(lemma_samples, k_clusters)
		return render_template('res_kmeans1.html', query=query,
		   x0_coordinates=x0_coordinates, y0_coordinates=y0_coordinates, z0_coordinates=z0_coordinates,
		   x1_coordinates=x1_coordinates, y1_coordinates=y1_coordinates, z1_coordinates=z1_coordinates,
		   x2_coordinates=x2_coordinates, y2_coordinates=y2_coordinates, z2_coordinates=z2_coordinates,
		   x3_coordinates=x3_coordinates, y3_coordinates=y3_coordinates, z3_coordinates=z3_coordinates,
		   x4_coordinates=x4_coordinates, y4_coordinates=y4_coordinates, z4_coordinates=z4_coordinates,
			titles0=titles0, titles1=titles1, titles2=titles2, titles3=titles3, titles4=titles4)
	else:
		return render_template('res_kmeans1.html', query=query)


#TODO: update for stacked barchart
@app.route('/res_stats/<query>', methods=["GET"]) #user statistics for iframe
def res_stats(query):
	pmid_list = query.split('+') #list of string pmids
	venn_data = make_venn(pmid_list)

	input_click_citations = statsSelfInfo(query)
	statistics = get_statistics(query)
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
	#app.run(host='0.0.0.0') if you are not running the app with uwsgi!


#TODO: Add unit tests and such for Git & Travis UI
########### GRAVEYARD ##########################################################


