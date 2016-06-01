from flask import Flask, render_template, request, flash, url_for, redirect, session, g
from flask_wtf import Form
from wtforms import StringField, TextField, SelectField
import sqlite3, gc, time, datetime
from content_management import runCosines
from cogeCrawled_db import connection
from Entrez_IR import getMainInfo, getCitationIDs, getCitedInfo, parsePMC, getContentPMC


#If running in a virtual enviornment, must have modules also (pip) installed in that virtualenv! 
#flask, Flask-WTF, nltk, bs4, lxml, requests, Bio



app = Flask(__name__)

#Content from content_management.py
cosines = runCosines() 



#Must return TWO things or else this error: #too many values to unpack (expected 2)
#Need to put everything in lists for output

#pmid "2342342" gave the error "list index out of range"

def run_IR_in_db(user_input):
	self_title, self_authors, self_journal = getMainInfo(user_input)
	pmc_ids = getCitationIDs(user_input)
	target_title, target_authors, target_journals, target_urls = getCitedInfo(pmc_ids)
	main_info = list(self_title)
	return main_info, target_title


def run_IR_not_db(user_input):
	self_title, self_authors, self_journal = getMainInfo(user_input)
	pmc_ids = getCitationIDs(user_input)
	target_title, target_authors, target_journals, target_urls = getCitedInfo(pmc_ids)
	main_info = list(self_title)
	getContentPMC(user_input, pmc_ids)
	return main_info, target_title



# def run_organismNER(user_prefix):
# 	latin_short, latin_long, wordnet_names = loadDocuments(user_prefix, 5)
# 	ner_info = list(zip(latin_short, latin_long, wordnet_names))
# 	return ner_info


#Main page
#User inputs a pubmed id and is then redirected to /results
#Prints sample results from coge publication
@app.route('/cogecrawl/', methods=["GET", "POST"])
def cogecrawl():
	error = None
	try:
		if request.method == "POST":
			attempted_pmid = request.form['pmid']
			#flash(attempted_pmid)
	except Exception as e:
		#flash(e)
		return render_template("dashboard.html", error=error) 
	return render_template('dashboard.html', cosines = cosines) #I should do the example page here :')





#Creat Form for handling user-entered pmid
#Need to pass pmid in form to MainCrawler.py
class pmidForm(Form):
	pmid = TextField('PubmedID')



#Getting Flask-WTFs to work with sqlite3 here
#This function uses user entry to run the MainCrawler.py 
#Prints the results under the "Information" tab
#User entered pmid is entered into sqlite3 database
#
#### TO DO: IF DONT RECOGNIZE PMID IN DB, THEN SCRAPE. IF RECOGNIZE, UPLOAD
#
@app.route('/results/', methods=["GET", "POST"])
def trying():
	form = pmidForm()
	try:
		if request.method == 'POST':
			pmid = form.pmid.data #THIS IS THE USER INPUT FROM THE FORM #referencing 'class pmidForm'
			user_input = str(pmid)


			#Connect to database
			conn, c = connection()
			#Check database for pmid
			c.execute("SELECT * FROM cogeCrawled WHERE pmids = (?)", (pmid, ))  #This code tells me 'db is locked'
			check1 = c.fetchone()
			#if the entry exists in the db already...
			

			#if the entry does NOT exist in the db already, will need to retireve text
			if check1 is None: 
				flash('congrats you entered a new pubmedid lol')
				#Using user_input for IR
				main_info, target_journals = run_IR_not_db(user_input)

				#add to sqlite3 database entry
				unix = time.time()
				date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
				conn, c = connection()
				c.execute("INSERT INTO cogeCrawled (datestamp, pmids) VALUES (?, ?)", (date, pmid)) #put user pmid into db
				conn.commit()
				flash("Writing PubmedID to database: successful")


			#if the entry IS in the db, no need to retireve text 
			if check1 is not None:
				flash("alreay exists in database :) ")
				#Do 
				main_info, target_journals = run_IR_in_db(user_input)
				#user_prefix = '/Users/hclent/Desktop/webdev-biotool/flask/'+user_input
				#ner_results = run_organismNER(user_prefix)
				#print(ner_results)


			#End cursor and connection to database
			c.close()
			conn.close()


			#Housekeeping
			gc.collect() #garbage collector for cleaning up unneeded stuff
			session['entered_id'] = True
			session['engaged'] = 'engaged' 
 

		return render_template('results.html', form=form, main_info = main_info, cosines=cosines)
	
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
# #Just runs pmc_spider and returns Journals for pmids in the cache
# def run_pmcSpider(user_input):
# 	titles, urls, authors = pmc_spider(1, user_input)
# 	flash("Got main info")
# 	main_info = list(zip(titles, authors, urls))
# 	return titles, urls, authors, main_info


# #Runs pmc_spider and get_texts for pmids not in the cache
# def run_FullMainCrawler(user_input):
# 	titles, urls, authors, main_info = run_pmcSpider(user_input)
# 	flash("Got main info... now texts...")
# 	journals = get_text(urls, user_input) #list
# 	flash("Got texts")
# 	return main_info, journals
