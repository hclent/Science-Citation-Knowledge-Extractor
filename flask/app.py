from flask import Flask, render_template, request, flash, url_for, redirect, session, g, Blueprint
from flask_wtf import Form
from wtforms import StringField, TextField, SelectField
import sqlite3, gc, time, datetime, pickle, os.path
import sys
from werkzeug.serving import run_simple
sys.path.append('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/')
from database_management import connection #mine
from content_management import * #mine
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
	error = None
	with open('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_mainfo.pickle', 'rb')as f:
		main_info = pickle.load(f)
	with open('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_journals.pickle', 'rb')as f:
		journals = pickle.load(f)
	return render_template("dashboard.html", main_info=main_info, journals=journals)



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
			print(pmid_list)

			q = '+'
			query = str(q.join(pmid_list))
			print("query: " + str(query))

			main_info = []
			target_journals = []
			target_dates = []
			data_samples = []
			ners = []

			for user_input in pmid_list:
				print(str(user_input))
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
					flash('congrats you entered a new pubmedid lol')
					#Using user_input for Information Retireval of "main info"
					self_info, main, journals, dates = run_IR_not_db(user_input)

					#add self_info to inputPapers database entry
					# for tup in self_info:
					# 	title = tup[0]
					# 	s = ', '
					# 	author = str(s.join(tup[1]))
					# 	journal = tup[2]
					# 	pubdate = tup[3]
					# 	url = tup[4]
                    #
					# 	unix = time.time()
					# 	date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
					# 	conn, c = connection()
					# 	c.execute("INSERT INTO inputPapers (datestamp, pmid, title, author, journal, pubdate, url) VALUES (?, ?, ?, ?, ?, ?, ?)", (date, user_input, title, author, journal, pubdate, url)) #put user pmid into db
					# 	conn.commit()
					# 	logging.info("Writing self_info to inputPapers db")
                    #
                    #
					# for tup in main:
					# 	logging.info("TUP IN MAIN: ")
					# 	logging.info(tup)
					# 	pmcid = tup[0]
					# 	title = tup[1]
					# 	s = ', '
					# 	author = str(s.join(tup[2]))
					# 	journal = tup[3]
					# 	pubdate = tup[4]
					# 	url = tup[5]
					# 	unix = time.time()
					# 	date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
					# 	conn, c = connection()
					# 	c.execute("INSERT INTO citations (datestamp, pmcid, title, author, journal, pubdate, citesPmid, url) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (date, pmcid, title, author, journal, pubdate, user_input, url)) #put user pmid into db
					# 	conn.commit()
					# 	logging.info("Writing main_info to citations db")

					for mi in main:
						main_info.append(mi)
					logging.info("done with main info list")
					for j in journals:
						target_journals.append(j)
					logging.info("done with journal list")
					lenjournals = (len(target_journals))
					logging.info("there are "+str(lenjournals)+" publications")
					for d in dates:
						target_dates.append(d)
					logging.info("done with dates list")


					logging.info("beginning multi-preprocessing")
					data, named_entities = do_ALL_multi_preprocessing(user_input)
					for d in data:
						data_samples.append(d)
					for n in named_entities:
						ners.append(n)

					## Once all the data has been acquired, do topic modeling
					## Only want to save final topic model (not running topic models)
					if user_input == pmid_list[-1]:

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


					#after successfully retrieving papers, annotating, and doing topic models,
					#add self_info to inputPapers database entry
					for tup in self_info:
						title = tup[0]
						s = ', '
						author = str(s.join(tup[1]))
						journal = tup[2]
						pubdate = tup[3]
						url = tup[4]

						unix = time.time()
						date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
						conn, c = connection()
						c.execute("INSERT INTO inputPapers (datestamp, pmid, title, author, journal, pubdate, url) VALUES (?, ?, ?, ?, ?, ?, ?)", (date, user_input, title, author, journal, pubdate, url)) #put user pmid into db
						conn.commit()
						logging.info("Writing self_info to inputPapers db")

					#add main to citations database table
					for tup in main:
						logging.info("TUP IN MAIN: ")
						logging.info(tup)
						pmcid = tup[0]
						title = tup[1]
						s = ', '
						author = str(s.join(tup[2]))
						journal = tup[3]
						pubdate = tup[4]
						url = tup[5]
						unix = time.time()
						date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
						conn, c = connection()
						c.execute("INSERT INTO citations (datestamp, pmcid, title, author, journal, pubdate, citesPmid, url) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (date, pmcid, title, author, journal, pubdate, user_input, url)) #put user pmid into db
						conn.commit()
						logging.info("Writing main_info to citations db")




				#if the entry IS in the db, no need to retrieve text from Entrez, just grab
				if check1 is not None:
					flash("alreay exists in database :) ")
					#Using user_input for Information Retireval of "main info"
					self_info, main, journals, dates = run_IR_in_db(user_input)

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


					data, named_entities = do_SOME_multi_preprocessing(user_input)
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
						journal_years = range_info[0]
						start_year = journal_years[0]
						end_year = journal_years[1]
						range_years = str(q.join(journal_years)) #2009+2016
						logging.info("range years: "+range_years)

						unique_publications = range_info[1]
						unique_journals = range_info[2]

						print(user_input+" is the last one (LSA)")
						print_lsa(query, user_input, jsonDict) #print lsa topic model to json
						print(user_input+" is the last one (LDA)")
						print_lda(query, user_input, jsonLDA) #print lda topic model to json


						print_data_and_nes(query, user_input, data_samples, ners) #print data_samples and nes_list to pickle


				#End cursor and connection to database
				c.close()
				conn.close()


				#Housekeeping
				gc.collect() #garbage collector for cleaning up unneeded stuff
				session['entered_id'] = True
				session['engaged'] = 'engaged'

		return render_template('results.html', form=form,
	   			main_info = main_info, target_journals = target_journals, query=query, range_years=range_years,
			   start_year=start_year, end_year=end_year, unique_publications=unique_publications, unique_journals=unique_journals)


	except Exception as e:
		return(str(e))




################ Forms for visualization toggle ################


class visOptions(Form):
	k_val = SelectField('k_val', choices=[(2,'k=2'),(3,'k=3'),(4,'k=4'),(5,'k=5')])
	w_words = SelectField('w_words', choices=[(4, 'w=4'),(5, 'w=5'),(6, 'w=6'), (7, 'w=7') ])


class nesOptions(Form):
	w_words = SelectField('w_words', choices=[(2, 'N'),(3, '3'),(10, '10'), (25, '25'), (50, '50'),(100, '100'),(200, '200'), (300, '300')])


################ Default CoGe Data #############################
@app.route('/cogelsa/', methods=["GET","POST"]) #default coge lsa for iframe
def cogelsa():
	form = visOptions(secret_key='super secret key')
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		print("the k value is " + str(k_clusters))
		query = '18952863+18269575'

		data_filename = "/home/hclent/data/18269575/data_samples_18952863+18269575.pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		print("rerunning the analysis")
		k = int(k_clusters)
		jsonLSA = run_lsa1(data_samples, k)
		print("did it all!")
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
		print("the k value is " + str(k_clusters))
		num_words = form.w_words.data
		print("the w value is "+str(num_words))
		query = '18952863+18269575'

		data_filename = "/home/hclent/data/18269575/data_samples_18952863+18269575.pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		print("rerunning the analysis")
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
	completeName = "/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/journal-publications.json"
	with open(completeName) as load_data:
		journals = json.load(load_data) #doesn't need to be parsed
	return render_template('coge_journals.html', journals=journals)


@app.route('/cogewordcloud/', methods=["GET","POST"]) #default coge NES Word Cloud for iframe
def cogewordcloud():
	form = nesOptions(secret_key='super secret key')
	if request.method == 'POST':

		nes_categories = request.form.getlist('check')
		print(nes_categories)
		w_number = form.w_words.data
		print("the w value is "+str(w_number))

		nes_list =  pickle.load(open("/home/hclent/data/18269575/nes_18952863+18269575.pickle", "rb")) #pre-processed already

		print(nes_list)
		wordcloud_data = vis_wordcloud(nes_list, nes_categories, w_number)
		return render_template('coge_wordcloud.html',  wordcloud_data=wordcloud_data)
	else:
		#Default data
		completeName = "/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_wcloud1.json"
		with open(completeName) as load_data:
			wordcloud_data = json.load(load_data) #this result doesn't need to be parsed but unsure how to write that in javascript
		#so i'm going to read it in as a string that needs to be parsed anyway
		wordcloud_data = re.sub('\'', '\"', str(wordcloud_data)) #json needs double quotes, not single quotes
		return render_template('coge_wordcloud.html', wordcloud_data=wordcloud_data)

#default data has only 1 of 2 pmids...
@app.route('/cogeheatmap/', methods=["GET","POST"]) #default coge NES heatmap for iframe
def cogeheatmap():
	form = nesOptions(secret_key='super secret key')
	if request.method == 'POST':
		nes_categories = request.form.getlist('check')
		print(nes_categories)
		w_number = form.w_words.data
		print("the w value is "+str(w_number))

		nes_list =  pickle.load(open("/home/hclent/data/18269575/nes_18952863+18269575.pickle", "rb")) #pre-processed already
		data_samples =  pickle.load(open("/home/hclent/data/18269575/data_samples_18952863+18269575.pickle", "rb")) #pre-processed already

		x_docs, y_words, z_counts = vis_heatmap(data_samples, nes_list, nes_categories, w_number)
		print(z_counts)
		print(x_docs)
		print(y_words)
		return render_template('coge_heatmap2.html', z_counts=z_counts, x_docs=x_docs, y_words=y_words)
	else:
		#Default data
		return render_template('coge_heatmap1.html')

#only using pmid 18269575
@app.route('/cogekmeans/', methods=["GET","POST"]) #default coge kmreans for iframe
def cogekmeans():
	form = visOptions(secret_key='super secret key')
	if request.method == 'POST':

		data_samples =  pickle.load(open("/home/hclent/data/18269575/data_samples_18952863+18269575.pickle", "rb")) #pre-processed already

		k_clusters = form.k_val.data #2,3,4,or 5
		print("the k value is " + str(k_clusters))
		x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates = vis_kmeans(data_samples, k_clusters)
		print(x0_coordinates)
		print(x1_coordinates)
		print(x2_coordinates)
		print(x3_coordinates)
		print(x4_coordinates)
		return render_template('coge_kmeans2.html', x0_coordinates=x0_coordinates, y0_coordinates=y0_coordinates, z0_coordinates=z0_coordinates,
							   x1_coordinates=x1_coordinates, y1_coordinates=y1_coordinates, z1_coordinates=z1_coordinates,
							   x2_coordinates=x2_coordinates, y2_coordinates=y2_coordinates, z2_coordinates=z2_coordinates,
							   x3_coordinates=x3_coordinates, y3_coordinates=y3_coordinates, z3_coordinates=z3_coordinates,
							   x4_coordinates=x4_coordinates, y4_coordinates=y4_coordinates, z4_coordinates=z4_coordinates)
	else:
		return render_template('coge_kmeans.html')


############### Results visualizations #########################################
@app.route('/resjournals/<query>/<range_years>', methods=["GET", "POST"]) #user journals for iframe
def resjournals(query, range_years):
	#need to get last user_input
	print("in routine res-journals")
	print("JOURNALS ID: " +str(query))
	print("YEARS RANGE: " +str(range_years))
	#Need years for range
	years_list = range_years.split('+')
	s_year = years_list[0]
	e_year = years_list[1]

	#only want to load the json for the LAST id in the query (so includes all)
	pmid_list = query.split('+') #list of string pmids
	last_entry = pmid_list[-1]
	print("the last entry is: " + str(last_entry))

	file_name = "journals_"+str(query)+".json"
	print("last entry's JOURNAL is named: " + str(file_name))
	savePath = "/home/hclent/data/"+str(last_entry)
	completeName = os.path.join(savePath, file_name)
	print("complete file: " + str(completeName))
	with open(completeName) as load_data:
		journals = json.load(load_data)
	print(journals) #str
	return render_template('results_journals.html', journals=journals, s_year=s_year, e_year=e_year)


@app.route('/reslsa/<query>', methods=["GET", "POST"]) #user lsa for iframe
def reslsa(query):
	form = visOptions(secret_key='super secret key')
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		print("the k value is " + str(k_clusters))
		# pmid_list = query.split('+') #list of string pmids
		# print(pmid_list)
		# data_samples = []
		# for pmid in pmid_list:
		# 	data, nes = do_SOME_multi_preprocessing(pmid)
		# 	for d in data:
		# 		data_samples.append(d)

		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		data_filename = "/home/hclent/data/"+str(last_entry)+"/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		num_pubs = int(len(data_samples))
		print("there are  "+str(num_pubs)+ " publications")
		print("rerunning the analysis")
		k = int(k_clusters)
		if num_pubs < k:
			logging.info("k value is larger than number of publications")
			#flash("For LSA, you cannot have more topics than documents. Try again")
		temp_jsonDict = run_lsa1(data_samples, k)
		print("did it all!")
		return render_template('results_lsa.html', query=query, jsonDict=temp_jsonDict)
	else:
		#need to get last user_input
		#use id to do stuff
		print("in routine resla")
		print("RES-LSA ID: " +str(query))
		#only want to load the json for the LAST id in the query (so includes all)
		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		print("the last entry is: " + str(last_entry))
		file_name = "lsa_"+str(query)+".json" #file named after query
		print("last entry's LSA is named: " + str(file_name))
		savePath = "/home/hclent/data/"+str(last_entry) #saved in folder of last pmid
		completeName = os.path.join(savePath, file_name)
		print("complete file: " + str(completeName))
		with open(completeName) as load_data:
			jsonDict = json.load(load_data)
		return render_template('results_lsa.html', query=query, jsonDict=jsonDict)


@app.route('/reslda/<query>', methods=["GET", "POST"]) #user lda for iframe
def reslda(query):
	form = visOptions(secret_key='super secret key')
	if request.method == 'POST':
		k_clusters = form.k_val.data #2,3,4,or 5
		print("the k value is " + str(k_clusters))
		num_words = form.w_words.data
		print("the w value is "+str(num_words))

		# data_samples = []
		# for pmid in pmid_list:
		# 	data, nes = do_SOME_multi_preprocessing(pmid)
		# 	for d in data:
		# 		data_samples.append(d)
		#

		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		data_filename = "/home/hclent/data/"+str(last_entry)+"/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		print("rerunning the analysis")
		k = int(k_clusters)
		w = int(num_words)
		temp_jsonLDA = run_lda1(data_samples, k, w)
		return render_template('coge_lda.html', form=form, jsonLDA=temp_jsonLDA, query=query)
	else:
		#need to get last user_input
		#use id to do stuff
		print("in routine reslDa")
		print("RES-LDA ID: " +str(query))
		#only want to load the json for the LAST id in the query (so includes all)
		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		print("the last entry is: " + str(last_entry))
		file_name = "lda_"+str(query)+".json" #file named after query
		print("last entry's LDA is named: " + str(file_name))
		savePath = "/home/hclent/data/"+str(last_entry) #saved in folder of last pmid
		completeName = os.path.join(savePath, file_name)
		print("complete file: " + str(completeName))
		with open(completeName) as load_data:
			jsonLDA = json.load(load_data)
		print(jsonLDA)
		return render_template('results_lda.html', form=form, jsonLDA=jsonLDA, query=query)


@app.route('/reswordcloud/<query>', methods=["GET", "POST"]) #user wordcloud for iframe
def reswordcloud(query):
	form = nesOptions(secret_key='super secret key')
	if request.method == 'POST':

		nes_categories = request.form.getlist('check')
		print(nes_categories)
		w_number = form.w_words.data
		print("the w value is "+str(w_number))

		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		filename = "/home/hclent/data/"+str(last_entry)+"/nes_"+str(query)+".pickle"
		nes_list =  pickle.load(open(filename, "rb")) #pre-processed already

		print(nes_list)
		wordcloud_data = vis_wordcloud(nes_list, nes_categories, w_number)
		return render_template('results_wordcloud.html', query=query,  wordcloud_data=wordcloud_data)
	else:
		nes_categories= ['BioProcess', 'CellLine', 'Cellular_component', 'Family', 'Gene_or_gene_product', 'Organ', 'Simple_chemical', 'Site', 'Species', 'TissueType']
		print(nes_categories)
		w_number = 10
		print("the w value is "+str(w_number))

		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		filename = "/home/hclent/data/"+str(last_entry)+"/nes_"+str(query)+".pickle"
		nes_list =  pickle.load(open(filename, "rb")) #pre-processed already

		wordcloud_data = vis_wordcloud(nes_list, nes_categories, w_number)

		return render_template('results_wordcloud.html', query=query, wordcloud_data=wordcloud_data)


@app.route('/res_heatmap/<query>', methods=["GET", "POST"]) #user heatmap for iframe
def res_heatmap(query):
	form = nesOptions(secret_key='super secret key')
	if request.method == 'POST':
		nes_categories = request.form.getlist('check')
		print(nes_categories)
		w_number = form.w_words.data
		print("the w value is "+str(w_number))

		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		nes_filename = "/home/hclent/data/"+str(last_entry)+"/nes_"+str(query)+".pickle"
		nes_list =  pickle.load(open(nes_filename, "rb")) #pre-processed already

		data_filename = "/home/hclent/data/"+str(last_entry)+"/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		x_docs, y_words, z_counts = vis_heatmap(data_samples, nes_list, nes_categories, w_number)
		print(z_counts)
		print(x_docs)
		print(y_words)
		return render_template('results_heatmap.html', query=query, z_counts=z_counts, x_docs=x_docs, y_words=y_words)
	else:
		nes_categories= ['BioProcess', 'CellLine', 'Cellular_component', 'Family', 'Gene_or_gene_product', 'Organ', 'Simple_chemical', 'Site', 'Species', 'TissueType']
		print(nes_categories)
		w_number = 10
		print("the w value is "+str(w_number))

		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		nes_filename = "/home/hclent/data/"+str(last_entry)+"/nes_"+str(query)+".pickle"
		nes_list =  pickle.load(open(nes_filename, "rb")) #pre-processed already

		data_filename = "/home/hclent/data/"+str(last_entry)+"/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		x_docs, y_words, z_counts = vis_heatmap(data_samples, nes_list, nes_categories, w_number)
		return render_template('results_heatmap.html', query=query, z_counts=z_counts, x_docs=x_docs, y_words=y_words)


@app.route('/res_kmeans/<query>', methods=["GET", "POST"]) #user k-means for iframe
def res_kmeans(query):
	form = visOptions(secret_key='super secret key')
	if request.method == 'POST':

		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		data_filename = "/home/hclent/data/"+str(last_entry)+"/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		k_clusters = form.k_val.data #2,3,4,or 5
		print("the k value is " + str(k_clusters))
		x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates = vis_kmeans(data_samples, k_clusters)
		print(x0_coordinates)
		print(x1_coordinates)
		print(x2_coordinates)
		print(x3_coordinates)
		print(x4_coordinates)
		return render_template('res_kmeans1.html', query=query,
		   x0_coordinates=x0_coordinates, y0_coordinates=y0_coordinates, z0_coordinates=z0_coordinates,
		   x1_coordinates=x1_coordinates, y1_coordinates=y1_coordinates, z1_coordinates=z1_coordinates,
		   x2_coordinates=x2_coordinates, y2_coordinates=y2_coordinates, z2_coordinates=z2_coordinates,
		   x3_coordinates=x3_coordinates, y3_coordinates=y3_coordinates, z3_coordinates=z3_coordinates,
		   x4_coordinates=x4_coordinates, y4_coordinates=y4_coordinates, z4_coordinates=z4_coordinates)
	else:
		pmid_list = query.split('+') #list of string pmids
		last_entry = pmid_list[-1]
		data_filename = "/home/hclent/data/"+str(last_entry)+"/data_samples_"+str(query)+".pickle"
		data_samples =  pickle.load(open(data_filename, "rb")) #pre-processed already

		k_clusters = 3 #default is 3
		print("the k value is " + str(k_clusters))
		x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, x2_coordinates, y2_coordinates, z2_coordinates, x3_coordinates, y3_coordinates, z3_coordinates, x4_coordinates, y4_coordinates, z4_coordinates = vis_kmeans(data_samples, k_clusters)
		print(x0_coordinates)
		print(x1_coordinates)
		print(x2_coordinates)
		print(x3_coordinates)
		print(x4_coordinates)
		return render_template('res_kmeans1.html', query=query,
		   x0_coordinates=x0_coordinates, y0_coordinates=y0_coordinates, z0_coordinates=z0_coordinates,
		   x1_coordinates=x1_coordinates, y1_coordinates=y1_coordinates, z1_coordinates=z1_coordinates,
		   x2_coordinates=x2_coordinates, y2_coordinates=y2_coordinates, z2_coordinates=z2_coordinates,
		   x3_coordinates=x3_coordinates, y3_coordinates=y3_coordinates, z3_coordinates=z3_coordinates,
		   x4_coordinates=x4_coordinates, y4_coordinates=y4_coordinates, z4_coordinates=z4_coordinates)

#################### OTHER ####################################################
@app.route('/testingstuff/')
def testingshit():
	print("testing stuff")
	return render_template('test.html')


#Handles 404 errors
@app.errorhandler(404)
def page_not_found(e):
	return("you shouldnt be here!!")


#Configuration settings
if __name__ == '__main__':
	run_simple('0.0.0.0', 5000, app, use_reloader=True)
	#app.run(host='0.0.0.0') #dont want app.run() for uwsgi

########### GRAVEYARD ##########################################################


