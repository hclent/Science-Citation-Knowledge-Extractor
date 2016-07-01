from flask import Flask, render_template, request, flash, url_for, redirect, session, g
from flask_wtf import Form
from wtforms import StringField, TextField, SelectField
import sqlite3, gc, time, datetime, pickle
from cogeCrawled_db import connection
from content_management import *
from processors import * 

#If running in a virtual enviornment, must have modules also (pip) installed in that virtualenv! 
#flask, Flask-WTF, nltk, bs4, lxml, requests, Biopython, pyProcessors


app = Flask(__name__)

# BUGS:
# Won't run code with only one citation
# Sometimes: "NCBI C++ Exception: Error: TXCLIENT(CException::eUnknown) ... Read failed: EOF (the other side has unexpectedly closed connection ...
	# 6-1-2016: 5 of these errors


#Main page
#User inputs a pubmed id and is then redirected to /results
#Prints sample results from 1 coge publication
@app.route('/cogecrawl/', methods=["GET", "POST"])
def cogecrawl():
	error = None
	with open('p18952863.pickle', 'rb')as f:
		main_info = pickle.load(f)
	with open('18952863_NES.p', 'rb')as f:
		ner = pickle.load(f)
	journals = ['BMC Genomics', 'Molecular Genetics and Genomics', 'Nature genetics', 'PLoS ONE', 'Frontiers in Plant Science', 'Frontiers in Plant Science', 'Frontiers in Genetics', 'PLoS ONE', 'BioMed Research International', 'Journal of Experimental Botany', 'Frontiers in Plant Science', 'Scientific Reports', 'Frontiers in Plant Science', 'PLoS ONE', 'BMC Genomics', 'GigaScience', 'Frontiers in Plant Science', 'PLoS ONE', 'Frontiers in Plant Science', 'BMC Evolutionary Biology', 'Nucleic Acids Research', 'Nucleic Acids Research', 'Frontiers in Plant Science', 'BMC Plant Biology', 'Plant Molecular Biology Reporter / Ispmb', 'BMC Genomics', 'Plant Physiology', 'BMC Genomics', 'BMC Genomics', 'The Plant journal : for cell and molecular biology', 'The Plant Cell', 'BMC Genomics', 'BMC Genomics', 'GigaScience', 'PLoS Genetics', 'Philosophical Transactions of the Royal Society B: Biological Sciences', 'Molecular Biology and Evolution', 'Frontiers in Plant Science', 'International Journal of Molecular Sciences', 'Scientific Reports', 'PLoS ONE', 'Journal of Experimental Botany', 'PLoS ONE', 'PLoS ONE', 'International Journal of Molecular Sciences', 'BMC Genetics', 'BMC Bioinformatics', 'BMC Bioinformatics', 'BMC Evolutionary Biology', 'BMC Genomics', 'BMC Genomics', 'Genome Biology and Evolution', 'BMC Genomics', 'Nucleic Acids Research', 'Genome Biology', 'Journal of Experimental Botany', 'PLoS ONE', 'PLoS ONE', 'Bioinformatics', 'Frontiers in Plant Science', 'BMC Bioinformatics', 'Frontiers in Plant Science', 'The Plant Cell', 'BMC Bioinformatics', 'Mobile Genetic Elements', 'Proceedings of the National Academy of Sciences of the United States of America', 'Frontiers in Plant Science', 'Frontiers in plant science', 'Nucleic Acids Research', 'PLoS ONE', 'Plant Physiology', 'The Plant Cell', 'Genome Biology and Evolution', 'Genome Biology', 'The Plant Cell', 'Nucleic Acids Research', 'BMC Plant Biology', 'PLoS ONE', 'BMC Evolutionary Biology', 'BMC Plant Biology', 'Annals of Botany', 'PLoS Biology', 'Journal of Molecular Evolution', 'Bioinformatics', 'Genome Biology and Evolution', 'Genome Research']
	try:
		if request.method == "POST":
			attempted_pmid = request.form['pmid']
			#flash(attempted_pmid)
	except Exception as e:
		#flash(e)
		return render_template("dashboard.html", error=error) 
	return render_template('dashboard.html', main_info=main_info, journals=journals, ner=ner) #I should do the example page here :')



#Creat Form for handling user-entered pmid
#Need to pass pmid in form to MainCrawler.py
class pmidForm(Form):
	pmid = TextField('PubmedID')



#Getting Flask-WTFs to work with sqlite3 here
#This function uses user entry to run the MainCrawler.py 
#Prints the results under the "Information" tab
#User entered pmid is entered into sqlite3 database
@app.route('/results/', methods=["GET", "POST"])
def trying():
	form = pmidForm()
	try:
		if request.method == 'POST':
			entry = form.pmid.data #THIS IS THE USER INPUT FROM THE FORM #referencing 'class pmidForm'
			pmid_list = multiple_pmid_input(entry)
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
				api = ProcessorsAPI(port=8886, keep_alive=True)

			
				#if the entry does NOT exist in the db already, will need to retireve text and annotate
				if check1 is None: 
					flash('congrats you entered a new pubmedid lol')
					#Using user_input for Information Retireval of "main info"
					main, journals = run_IR_not_db(user_input)
					for mi in main:
						main_info.append(mi)
					for j in journals:
						target_journals.append(j)	
					user_prefix = '/Users/hclent/Desktop/webdev-biotool/flask/data/'+user_input
					num = len(journals) #how many docs there are? 
					#Do_preprocessing() annotates the docs with pyProcessors  #Function in content_management.py
					data, named_entities = do_preprocessing(num, user_input, api)
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
					num = len(journals) #how many docs there are
					user_prefix = '/Users/hclent/Desktop/webdev-biotool/flask/data/'+user_input
					#Do_preprocessing() annotates the docs with pyProcessors 
					data, named_entities = already_have_preproc(num, user_input)
					for d in data:
						data_samples.append(d)
					for n in named_entities:
						ners.append(n)
					#Do visualization
					jsonDict = run_lsa1(user_input, data_samples, 2)


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



#Built in Flask debugger
if __name__ == '__main__':
	app.secret_key = 'super secret key'
	app.run(debug=True)


########### GRAVEYARD ########

# @app.route('/visdev/') #this is where I'm experimenting with data visualization
# def visDEV():
# 	return render_template('vis.html') 

