from flask import Flask, render_template, request, flash, url_for, redirect, session, g, Blueprint
from flask_wtf import Form
from wtforms import TextField, SelectField
import sqlite3, gc, time, datetime, pickle, os.path
import sys
from werkzeug.serving import run_simple
sys.path.append('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/')
from database_management import connection #mine
from content_management import * #mine
from citation_venn import make_venn #mine
from processors import *

#If running in a virtual enviornment, must have modules also (pip) installed in that virtualenv! 
#Flask-WTF, Biopython==1.67, py-Processors
#flask, nltk, bs4, lxml, requests

app = Flask(__name__, static_url_path='/hclent/Webdev-for-bioNLP-lit-tool/flask/static')
app.debug = True
app.secret_key = 'super secret key'


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
	form = pmidForm(secret_key='super secret key')
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

			main_info = [] #main_info are formatted citations from citations table in db
			target_journals = []
			target_dates = []
			target_urls = []
			data_samples = []
			ners = []


			for user_input in pmid_list:
				logging.info(str(user_input))
				user_input = str(user_input)

				############################################
				#Connect to database
				conn, c = connection()
				#Check database for pmid #Does the entry exists in the db already?
				#Check inputPapers instead of cogeCrawled
				c.execute("SELECT * FROM inputPapers WHERE pmid = (?)", (user_input, ))

				check1 = c.fetchone()


				#if the entry does NOT exist in the db already, will need to retrieve text and annotate
				if check1 is None:
					flash('new pubmedid lol')
					#Using user_input for Information Retireval of "main info"
					self_info, new_info, journals, dates, num_citations = run_IR_not_db(user_input)

					'''
					'main' renamed to 'new_info'
					need 'new_info' written to the database
					but need this information formatted to citation form before printing it as `main`
					so appending to 'main' will be moved to after new_info is put into db



					# for mi in main:
					# 	main_info.append(mi)
					# logging.info("done with main info list")
					'''
					for j in journals:
						target_journals.append(j)
					logging.info("done with journal list")
					lenjournals = (len(target_journals))
					logging.info("there are "+str(lenjournals)+" publications")
					for d in dates:
						target_dates.append(d)
					logging.info("done with dates list")

					# add main/new to citations database table
					# need these added to the db in order to annotate the papers

					for tup in new_info:
						logging.info("TUP IN MAIN: ")
						logging.info(tup)
						pmcid = tup[0]
						title = tup[1]
						s = ', '
						author = str(s.join(tup[2]))
						journal = tup[3]
						pubdate = tup[4]
						url = tup[5]
						abstract = tup[6]
						whole = tup[7]

						#can no longer write sents and tokens here because new_info must be added to db before annotation
						#sents = total_sentences[i]
						#tokens = sum_tokens[i]

						unix = time.time()
						date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
						conn, c = connection()
						c.execute(
							"INSERT INTO citations (datestamp, pmcid, title, author, journal, pubdate, citesPmid, url, abstract, whole_article) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
							(date, pmcid, title, author, journal, pubdate, user_input, url, abstract, whole))  # put user pmid into db
						conn.commit()

					logging.info("Writing citing new_info to citations db")

					logging.info("beginning multi-preprocessing")
					data, named_entities, total_sentences, sum_tokens = do_ALL_multi_preprocessing(user_input)

					#put total_sents and sum_tokens info into db
					i = 0
					for tup in new_info:
						logging.info(tup)
						pmcid = tup[0]
						logging.info(pmcid)
						sents = total_sentences[i]
						logging.info(sents)
						tokens = sum_tokens[i]
						logging.info(tokens)
						conn, c = connection()
						c.execute("UPDATE citations SET sents=?, tokens=? WHERE pmcid=?", (sents, tokens, pmcid))
						conn.commit()
						i += 1
					logging.info("updated the db")

					for d in data:
						data_samples.append(d)
					for n in named_entities:
						ners.append(n)


					## Once all the data has been acquired, (no topic modeling yet)
					## Only want to save final topic model (not running topic models)

					for tup in self_info:
						title = tup[0]
						s = ', '
						author = str(s.join(tup[1]))
						journal = tup[2]
						pubdate = tup[3]
						url = tup[4]
						#needs "num_citations"
						unix = time.time()
						date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
						conn, c = connection()
						c.execute("INSERT INTO inputPapers (datestamp, pmid, title, author, journal, pubdate, url, num_citations) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (date, user_input, title, author, journal, pubdate, url, num_citations)) #put user pmid into db
						conn.commit()
					logging.info("Writing self_info to inputPapers db")

					# If its NOT in the db, we also need to try scraping it and writing that info to db
					scrape_and_write_Input(user_input)


					if user_input == pmid_list[-1]: #if its the last pmid

						#Do Latent Semantic Analysis and return jsonDict for data vis
						jsonDict = run_lsa1(data_samples, 2)

						#Do Latent Dirichlet Allocation
						jsonLDA = run_lda1(data_samples, 3, 5)

						logging.info(user_input+" is the last one (JOURNALS)")
						range_info = print_journalvis(target_journals, target_dates, user_input, query) #e.g. [(2008, 2009), 10, 7]
						journal_years = range_info[0]
						start_year = journal_years[0]
						end_year = journal_years[1]
						range_years = str(q.join(journal_years))
						logging.info("range years: "+range_years)


						unique_publications = range_info[1]
						unique_journals = range_info[2]

						logging.info(user_input+" is the last one (LSA)")
						print_lsa(query, user_input, jsonDict) #print lsa topic model to json
						logging.info(user_input+" is the last one (LDA)")
						print_lda(query, user_input, jsonLDA) #print lda topic model to json


						print_data_and_nes(query, user_input, data_samples, ners) #print data_samples and nes_list to pickle


					#after info written to db, now can access db and get formated main_info (main)
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
					#Using user_input for Information Retireval of "main info"
					self_info, main, journals, dates, db_urls = run_IR_in_db(user_input)

					for mi in main:
						main_info.append(mi)
					logging.info("done with main info list")
					for j in journals:
						target_journals.append(j)
					logging.info("done with journals list")
					lenjournals = (len(target_journals))
					logging.info("there are "+str(lenjournals)+" publications")
					for d in dates:
						target_dates.append(d)
					logging.info("done with dates list")
					for url in db_urls:
						target_urls.append(url)
					logging.info("done with url list")


					data, named_entities, total_sentences, sum_tokens = do_SOME_multi_preprocessing(user_input)
					for d in data:
						data_samples.append(d)
					for n in named_entities:
						ners.append(n)


					#### errr, if its in the Db, its probably already been topic modeled. so instead should try to load the data,
					## and if can't load the data, do the analysis

					## Now that we have all the data, do the topic model
					## Only want to save final topic model (not running topic model)
					if user_input == pmid_list[-1]:

						# #Do visualization and Topic Modeling
						jsonDict = run_lsa1(data_samples, 2) #default = 2 "topics" for right now

						#latent dirichlet allocation
						jsonLDA = run_lda1(data_samples, 3, 5)

						logging.info(user_input+" is the last one (JOURNALS)")
						range_info = print_journalvis(target_journals, target_dates, user_input, query)
						logging.info(range_info)
						journal_years = range_info[0]
						start_year = journal_years[0]
						end_year = journal_years[1]
						range_years = str(q.join(journal_years)) #2009+2016
						logging.info("range years: "+range_years)

						unique_publications = range_info[1]
						unique_journals = range_info[2]

						logging.info(user_input+" is the last one (LSA)")
						print_lsa(query, user_input, jsonDict) #print lsa topic model to json
						logging.info(user_input+" is the last one (LDA)")
						print_lda(query, user_input, jsonLDA) #print lda topic model to json

						print_data_and_nes(query, user_input, data_samples, ners) #print data_samples and nes_list to pickle


				#End cursor and connection to database
				c.close()
				conn.close()


				#Housekeeping
				gc.collect() #garbage collector for cleaning up unneeded stuff
				session['entered_id'] = True
				session['engaged'] = 'engaged'


		citations_with_links = list(zip(main_info, target_urls))


		return render_template('results.html', form=form, citations_with_links=citations_with_links,
	   			main_info = main_info, target_journals = target_journals, query=query, range_years=range_years,
			   start_year=start_year, end_year=end_year, unique_publications=unique_publications, unique_journals=unique_journals)


	except Exception as e:
		return(str(e))




################ Forms for visualization toggle ################
class visOptions(Form):
	k_val = SelectField('k_val', choices=[(2,'k=2'),(3,'k=3'),(4,'k=4'),(5,'k=5')])
	w_words = SelectField('w_words', choices=[(4, 'w=4'),(5, 'w=5'),(6, 'w=6'), (7, 'w=7') ])


class nesOptions(Form):
	w_words = SelectField('w_words', choices=[(2, 'N'),(3, '3'),(10, '10'), (25, '25'), (50, '50'),(100, '100'),(200, '200'),
											  (300, '300')])

class corpusOptions(Form):
	corpus = SelectField('corpus', choices=[('startrek', 'startrek'),('frankenstein', 'frankenstein'),('youth', 'youth'),
											('darwin', 'darwin'), ('austen', 'austen'), ('paper1', 'paper1'),
											('paper2', 'paper2'),('paper3', 'paper3'), ('paper4', 'paper4'), ('paper5', 'paper5')])

################ Default CoGe Data #############################
@app.route('/cogelsa/', methods=["GET","POST"]) #default coge lsa for iframe
def cogelsa():
	form = visOptions(secret_key='super secret key')
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		query = '18952863+18269575'

		data_filename = "/home/hclent/data/data_samples/data_samples_18952863+18269575.pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		logging.info("rerunning the analysis")
		k = int(k_clusters)
		jsonLSA = run_lsa1(data_samples, k)
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
	form = visOptions(secret_key='super secret key')
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		num_words = form.w_words.data
		logging.info("the w value is "+str(num_words))
		query = '18952863+18269575'

		data_filename = "/home/hclent/data/data_samples/data_samples_18952863+18269575.pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		logging.info("rerunning the analysis")
		k = int(k_clusters)
		w = int(num_words)
		jsonLDA = run_lda1(data_samples, k, w)
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
	form = nesOptions(secret_key='super secret key')
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
	form = nesOptions(secret_key='super secret key')
	if request.method == 'POST':
		nes_categories = request.form.getlist('check')
		logging.info(nes_categories)
		w_number = form.w_words.data
		logging.info("the w value is "+str(w_number))

		nes_list =  pickle.load(open("/home/hclent/data/nes/nes_18952863+18269575.pickle", "rb")) #pre-processed already
		data_samples =  pickle.load(open("/home/hclent/data/data_samples/data_samples_18952863+18269575.pickle", "rb")) #pre-processed already

		x_docs, y_words, z_counts = vis_heatmap(data_samples, nes_list, nes_categories, w_number)
		#print(z_counts)
		#print(x_docs)
		#print(y_words)
		return render_template('coge_heatmap2.html', z_counts=z_counts, x_docs=x_docs, y_words=y_words)
	else:
		#Default data
		return render_template('coge_heatmap1.html')



@app.route('/cogekmeans/', methods=["GET","POST"]) #default coge k-means clustering for iframe
def cogekmeans():
	form = visOptions(secret_key='super secret key')
	if request.method == 'POST':

		data_samples =  pickle.load(open("/home/hclent/data/data_samples/data_samples_18952863+18269575.pickle", "rb")) #pre-processed already
		pmid_list = ['18952863', '18269575']
		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates, titles0, titles1, titles2, titles3, titles4 = vis_kmeans(data_samples, k_clusters, pmid_list)

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
	form = corpusOptions(secret_key='super secret key')
	#decide eligible_papers
	eligible_papers = [('paper1', '18952863', '/home/hclent/data/pmcids/259/367/2593677.txt')]
	if request.method == 'POST':
		logging.info("posted a thing in scifi!")
		corpus = form.corpus.data
		logging.info(corpus)
		query = "18952863+18269575"
		if corpus == 'startrek':
			color = 'rgb(63, 100, 168)'
			title = 'Star Trek: The Next Generation'
		if corpus == 'frankenstein':
			color = 'rgb(92, 59, 107)'
			title = 'Frankenstein; or, The Modern Prometheus'
		if corpus == 'youth':
			color = 'rgb(142, 7, 7)'
			title = 'Youth by Isaac Asimov'
		if corpus == 'darwin':
			color = 'rgb(8, 114, 32)'
			title = 'On The Origin of Species'
		if corpus == 'austen':
			color = 'rgb(191, 110, 167)'
			title = 'Pride and Prejudice'
		if corpus == 'paper1':
			color =  'rgb(8, 114, 32)'
			title =  'PMID: 18952863'
		# if corpus == 'paper2':
		# 	color =  'rgb(8, 114, 32)'
		# 	title = 'PMID: 18269575'
		x, y, names = vis_scifi(corpus, query, eligible_papers)

		return render_template('coge_scifi2.html', x=x, y=y, title=title, color=color, names=names, eligible_papers=eligible_papers)
	else:
		flash('Some input paper(s) are not avaliable')
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


@app.route('/reslsa/<query>', methods=["GET", "POST"]) #user lsa for iframe
def reslsa(query):
	form = visOptions(secret_key='super secret key')
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
	form = visOptions(secret_key='super secret key')
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
		logging.info("complete file: " + str(completeName))
		with open(completeName) as load_data:
			jsonLDA = json.load(load_data)
		logging.info(jsonLDA)
		return render_template('results_lda.html', form=form, jsonLDA=jsonLDA, query=query)


@app.route('/reswordcloud/<query>', methods=["GET", "POST"]) #user wordcloud for iframe
def reswordcloud(query):
	form = nesOptions(secret_key='super secret key')
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
	form = nesOptions(secret_key='super secret key')
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
		# print(z_counts)
		# print(x_docs)
		# print(y_words)
		popup = ' '
		return render_template('results_heatmap.html', query=query, z_counts=z_counts, x_docs=x_docs, y_words=y_words, popup=popup)
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

		popup = '<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button><strong>[ ! ]</strong> Default: N=10, from all categories.</div>'
		return render_template('results_heatmap.html', query=query, z_counts=z_counts, x_docs=x_docs, y_words=y_words, popup=popup)


@app.route('/res_kmeans/<query>', methods=["GET", "POST"]) #user k-means for iframe
def res_kmeans(query):
	form = visOptions(secret_key='super secret key')
	if request.method == 'POST':

		pmid_list = query.split('+') #list of string pmids
		data_filename = "/home/hclent/data/data_samples/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		k_clusters = form.k_val.data #2,3,4,or 5
		logging.info("the k value is " + str(k_clusters))
		x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates, titles0, titles1, titles2, titles3, titles4 = vis_kmeans(data_samples, k_clusters, pmid_list)
		# print(x0_coordinates)
		# print(x1_coordinates)
		# print(x2_coordinates)
		# print(x3_coordinates)
		# print(x4_coordinates)
		return render_template('res_kmeans1.html', query=query,
		   x0_coordinates=x0_coordinates, y0_coordinates=y0_coordinates, z0_coordinates=z0_coordinates,
		   x1_coordinates=x1_coordinates, y1_coordinates=y1_coordinates, z1_coordinates=z1_coordinates,
		   x2_coordinates=x2_coordinates, y2_coordinates=y2_coordinates, z2_coordinates=z2_coordinates,
		   x3_coordinates=x3_coordinates, y3_coordinates=y3_coordinates, z3_coordinates=z3_coordinates,
		   x4_coordinates=x4_coordinates, y4_coordinates=y4_coordinates, z4_coordinates=z4_coordinates,
			titles0=titles0, titles1=titles1, titles2=titles2, titles3=titles3, titles4=titles4)
	else:
		#pmid_list = query.split('+') #list of string pmids
		#last_entry = pmid_list[-1]
		#data_filename = "/home/hclent/data/"+str(last_entry)+"/data_samples_"+str(query)+".pickle"
		#data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		# k_clusters = 3 #default is 3
		# print("the k value is " + str(k_clusters))
		# x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates = vis_kmeans(data_samples, k_clusters)
		# print(x0_coordinates)
		# print(x1_coordinates)
		# print(x2_coordinates)
		# print(x3_coordinates)
		# print(x4_coordinates)
		# return render_template('res_kmeans1.html', query=query,
		#    x0_coordinates=x0_coordinates, y0_coordinates=y0_coordinates, z0_coordinates=z0_coordinates,
		#    x1_coordinates=x1_coordinates, y1_coordinates=y1_coordinates, z1_coordinates=z1_coordinates,
		#    x2_coordinates=x2_coordinates, y2_coordinates=y2_coordinates, z2_coordinates=z2_coordinates,
		#    x3_coordinates=x3_coordinates, y3_coordinates=y3_coordinates, z3_coordinates=z3_coordinates,
		#    x4_coordinates=x4_coordinates, y4_coordinates=y4_coordinates, z4_coordinates=z4_coordinates)
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
	form = corpusOptions(secret_key='super secret key')
	pmid_list = query.split('+')  # list of string pmids
	#decide eligible papers:
	eligible_papers = inputEligible(query)
	logging.info("eligible papers: " +str(eligible_papers))
	if request.method == 'POST':
		logging.info("posted a thing in scifi!")
		corpus = form.corpus.data
		logging.info(corpus)
		if corpus == 'startrek':
			color = 'rgb(63, 100, 168)'
			title = 'Star Trek: The Next Generation'
		if corpus == 'frankenstein':
			color = 'rgb(92, 59, 107)'
			title = 'Frankenstein; or, The Modern Prometheus'
		if corpus == 'youth':
			color = 'rgb(142, 7, 7)'
			title = 'Youth by Isaac Asimov'
		if corpus == 'darwin':
			color = 'rgb(8, 114, 32)'
			title = 'On The Origin of Species'
		if corpus == 'austen':
			color = 'rgb(191, 110, 167)'
			title = 'Pride and Prejudice'
		if corpus == 'paper1':
			color =  'rgb(8, 114, 32)'
			title =  str('PMID: '+ str(eligible_papers[0][1]))
		if corpus == 'paper2':
			color =  'rgb(8, 114, 32)'
			title = str('PMID: '+ str(eligible_papers[1][1]))
		if corpus == 'paper3':
			color = 'rgb(8, 114, 32)'
			title = str('PMID: ' + str(eligible_papers[2][1]))
		if corpus == 'paper4':
			color = 'rgb(8, 114, 32)'
			title = str('PMID: ' + str(eligible_papers[3][1]))
		if corpus == 'paper5':
			color = 'rgb(8, 114, 32)'
			title = str('PMID: ' + str(eligible_papers[4][1]))
		x, y, names = vis_scifi(corpus, query, eligible_papers)
		return render_template('results_scifi.html', x=x, y=y, title=title, color=color, query=query, names=names, eligible_papers=eligible_papers)
	else:
		logging.info("scifi analysis")
		corpus = 'startrek'
		title = 'Star Trek: The Next Generation'
		color = 'rgb(63, 100, 168)'
		x, y, names = vis_scifi(corpus, query, eligible_papers)
		logging.info("done with x and y")
		if len(eligible_papers) < len(pmid_list):
			flash('Some input paper(s) are not avaliable')
		return render_template('results_scifi.html', x=x, y=y, title=title, color=color, query=query, names=names, eligible_papers=eligible_papers)

#################### OTHER ####################################################
@app.route('/testingstuff/')
def testingshit():
	#print("testing stuff")
	completeName = "/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/templates/data.json"
	with open(completeName) as load_data:
		data = json.load(load_data) #doesn't need to be parsed but unsure how to write that in javascript
	return render_template('test.html', data=data)


#Handles 404 errors
@app.errorhandler(404)
def page_not_found(e):
	return("you shouldnt be here!!")


#Configuration settings
if __name__ == '__main__':
	run_simple('0.0.0.0', 5000, app, use_reloader=True)
	#app.run(host='0.0.0.0') #dont want app.run() for uwsgi

########### GRAVEYARD ##########################################################


