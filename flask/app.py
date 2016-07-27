from flask import Flask, render_template, request, flash, url_for, redirect, session, g, Blueprint
from flask_wtf import Form
from wtforms import StringField, TextField, SelectField
import sqlite3, gc, time, datetime, pickle
import sys
from werkzeug.serving import run_simple
sys.path.append('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/')
from cogeCrawled_db import connection
from content_management import *
from processors import * 

#If running in a virtual enviornment, must have modules also (pip) installed in that virtualenv! 
#flask, Flask-WTF, nltk, bs4, lxml, requests, Biopython, pyProcessors


app = Flask(__name__, static_url_path='/hclent/Webdev-for-bioNLP-lit-tool/flask/static')
app.debug = True
app.secret_key = 'super secret key'




#Main page
#Prints sample results from 2 coge publications
#User inputs a pubmed id and is then redirected to /results
@app.route("/cogecrawl/", methods=["GET", "POST"])
def cogecrawl():
	error = None
	with open('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_mainfo.pickle', 'rb')as f:
		main_info = pickle.load(f)
	with open('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_nes.pickle', 'rb')as f:
		ner = pickle.load(f)
	with open('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_journals.pickle', 'rb')as f:
		journals = pickle.load(f)
	try:
		if request.method == "POST":
			attempted_pmid = request.form['pmid']
			#flash(attempted_pmid)
	except Exception as e:
		#flash(e)
		return render_template("dashboard.html", error=error) 
	return render_template("dashboard.html", main_info=main_info, journals=journals, ner=ner)



#Creat Form for handling user-entered pmid
#Need to pass pmid in form to Entrez_IR.py
class pmidForm(Form):
	pmid = TextField('PubmedID')



#Getting Flask-WTFs to work with sqlite3 here
#This function uses user entry to run the Entrez_IR.py 
#User entered pmid is entered into sqlite3 database
@app.route("/results/", methods=["GET", "POST"])
def trying():
	form = pmidForm(secret_key='super secret key')
	try:
		if request.method == 'POST':
			entry = form.pmid.data #THIS IS THE USER INPUT FROM THE FORM #referencing 'class pmidForm'
			pmid_list = multiple_pmid_input(entry) #list for handling multiple pmids
			print(pmid_list)

			main_info = []
			target_journals = []
			data_samples = []
			ners = []

			for user_input in pmid_list:
				print(str(user_input))
				user_input = str(user_input)

				############################################
				#Connect to database
				conn, c = connection()
				#Check database for pmid #Does the entry exists in the db already?
				c.execute("SELECT * FROM cogeCrawled WHERE pmids = (?)", (user_input, ))
				check1 = c.fetchone()
				
				
				#Connect to Processors service
				api = connect_to_Processors(4343) #method from content_management


				#if the entry does NOT exist in the db already, will need to retireve text and annotate
				if check1 is None: 
					flash('congrats you entered a new pubmedid lol')
					#Using user_input for Information Retireval of "main info"
					main, journals = run_IR_not_db(user_input)
					for mi in main:
						main_info.append(mi)
					for j in journals:
						target_journals.append(j)

					data, named_entities = do_ALL_multi_preprocessing(user_input, api)
					for d in data:
						data_samples.append(d)
					for n in named_entities:
						ners.append(n)
					#Do Latent Semantic Analysis and return jsonDict for data vis
					jsonDict = run_lsa1(user_input, data_samples, 2)

					#add to sqlite3 database entry
					unix = time.time()
					date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
					conn, c = connection()
					c.execute("INSERT INTO cogeCrawled (datestamp, pmids) VALUES (?, ?)", (date, user_input)) #put user pmid into db
					conn.commit()
					flash("Writing PubmedID to database: successful")
					
				#if the entry IS in the db, no need to retireve text from Entrez, just grab  
				if check1 is not None:
					flash("alreay exists in database :) ")
					#Using user_input for Information Retireval of "main info"
					main, journals = run_IR_in_db(user_input)
					for mi in main:
						main_info.append(mi)
					for j in journals:
						target_journals.append(j)
						
					data, named_entities = do_SOME_multi_preprocessing(user_input, api)
					for d in data:
						data_samples.append(d)
					for n in named_entities:
						ners.append(n)
					# #Do visualization and Topic Modeling
					jsonDict = run_lsa1(user_input, data_samples, 2) #default = 2 "topics" for right now


				#End cursor and connection to database
				c.close()
				conn.close()


				#Housekeeping
				gc.collect() #garbage collector for cleaning up unneeded stuff
				session['entered_id'] = True
				session['engaged'] = 'engaged' 
	 

		return render_template('results.html', form=form, main_info = main_info, target_journals = target_journals, ners=ners, jsonDict=jsonDict)
	except Exception as e:
		return(str(e))



#Handles 404 errors
@app.errorhandler(404)
def page_not_found(e):
	return("you shouldnt be here!!")


#Configuration settings
if __name__ == '__main__':
	run_simple('0.0.0.0', 5000, app, use_reloader=True)
	#app.run(host='0.0.0.0') #dont want app.run() for uwsgi

########### GRAVEYARD ########

@app.route('/visdev/') #this is where I'm experimenting with data visualization
def visDEV():
	return render_template('vis_lsa.html')

