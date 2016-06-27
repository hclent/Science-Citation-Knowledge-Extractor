from flask import Flask, render_template, request, flash, url_for, redirect, session, g
from flask_wtf import Form
from wtforms import StringField, TextField, SelectField
import sqlite3, gc, time, datetime, pickle
from cogeCrawled_db import connection
from Entrez_IR import getMainInfo, getCitationIDs, getCitedInfo, parsePMC, getContentPMC
from organismNER import loadDocuments

#If running in a virtual enviornment, must have modules also (pip) installed in that virtualenv! 
#flask, Flask-WTF, nltk, bs4, lxml, requests, Biopython


app = Flask(__name__)


# BUGS:
# Won't run code with only one citation
# Sometimes: "NCBI C++ Exception: Error: TXCLIENT(CException::eUnknown) ... Read failed: EOF (the other side has unexpectedly closed connection ...
	# 6-1-2016: 5 of these errors


#Should switch these to content manager.... 
def run_IR_in_db(user_input):
	self_info = getMainInfo(user_input)
	pmc_ids = getCitationIDs(user_input)
	target_title, target_authors, target_journals, target_urls = getCitedInfo(pmc_ids)
	main_info = list(zip(target_title, target_authors, target_urls))
	return main_info, target_journals


def run_IR_not_db(user_input): #Not in db? Scrapes documents 
	self_info = getMainInfo(user_input)
	pmc_ids = getCitationIDs(user_input)
	target_title, target_authors, target_journals, target_urls = getCitedInfo(pmc_ids)
	main_info = list(zip(target_title, target_authors, target_urls))
	#Get XML 
	getContentPMC(user_input, pmc_ids)
	return main_info, target_journals



@app.route('/visdev/') #this is where I'm experimenting with data visualization
def visDEV():
	return render_template('vis.html') 

#Main page
#User inputs a pubmed id and is then redirected to /results
#Prints sample results from coge publication
@app.route('/cogecrawl/', methods=["GET", "POST"])
def cogecrawl():
	error = None
	with open('p18952863.pickle', 'rb')as f:
		main_info = pickle.load(f)
	journals = ['BMC Genomics', 'Molecular Genetics and Genomics', 'Nature genetics', 'PLoS ONE', 'Frontiers in Plant Science', 'Frontiers in Plant Science', 'Frontiers in Genetics', 'PLoS ONE', 'BioMed Research International', 'Journal of Experimental Botany', 'Frontiers in Plant Science', 'Scientific Reports', 'Frontiers in Plant Science', 'PLoS ONE', 'BMC Genomics', 'GigaScience', 'Frontiers in Plant Science', 'PLoS ONE', 'Frontiers in Plant Science', 'BMC Evolutionary Biology', 'Nucleic Acids Research', 'Nucleic Acids Research', 'Frontiers in Plant Science', 'BMC Plant Biology', 'Plant Molecular Biology Reporter / Ispmb', 'BMC Genomics', 'Plant Physiology', 'BMC Genomics', 'BMC Genomics', 'The Plant journal : for cell and molecular biology', 'The Plant Cell', 'BMC Genomics', 'BMC Genomics', 'GigaScience', 'PLoS Genetics', 'Philosophical Transactions of the Royal Society B: Biological Sciences', 'Molecular Biology and Evolution', 'Frontiers in Plant Science', 'International Journal of Molecular Sciences', 'Scientific Reports', 'PLoS ONE', 'Journal of Experimental Botany', 'PLoS ONE', 'PLoS ONE', 'International Journal of Molecular Sciences', 'BMC Genetics', 'BMC Bioinformatics', 'BMC Bioinformatics', 'BMC Evolutionary Biology', 'BMC Genomics', 'BMC Genomics', 'Genome Biology and Evolution', 'BMC Genomics', 'Nucleic Acids Research', 'Genome Biology', 'Journal of Experimental Botany', 'PLoS ONE', 'PLoS ONE', 'Bioinformatics', 'Frontiers in Plant Science', 'BMC Bioinformatics', 'Frontiers in Plant Science', 'The Plant Cell', 'BMC Bioinformatics', 'Mobile Genetic Elements', 'Proceedings of the National Academy of Sciences of the United States of America', 'Frontiers in Plant Science', 'Frontiers in plant science', 'Nucleic Acids Research', 'PLoS ONE', 'Plant Physiology', 'The Plant Cell', 'Genome Biology and Evolution', 'Genome Biology', 'The Plant Cell', 'Nucleic Acids Research', 'BMC Plant Biology', 'PLoS ONE', 'BMC Evolutionary Biology', 'BMC Plant Biology', 'Annals of Botany', 'PLoS Biology', 'Journal of Molecular Evolution', 'Bioinformatics', 'Genome Biology and Evolution', 'Genome Research']
	try:
		if request.method == "POST":
			attempted_pmid = request.form['pmid']
			#flash(attempted_pmid)

	except Exception as e:
		#flash(e)
		return render_template("dashboard.html", error=error) 
	return render_template('dashboard.html', main_info=main_info, journals=journals) #I should do the example page here :')



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
				num = len(target_journals) #how many docs there are

				user_prefix = '/Users/hclent/Desktop/webdev-biotool/flask/'+user_input


				#add to sqlite3 database entry
				unix = time.time()
				date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
				conn, c = connection()
				c.execute("INSERT INTO cogeCrawled (datestamp, pmids) VALUES (?, ?)", (date, pmid)) #put user pmid into db
				conn.commit()
				flash("Writing PubmedID to database: successful")


			#if the entry IS in the db, no need to retireve text from Entrez, just grab  
			if check1 is not None:
				flash("alreay exists in database :) ")
				#Do  
				main_info, target_journals = run_IR_in_db(user_input)
				num = len(target_journals) #how many docs there are
				user_prefix = '/Users/hclent/Desktop/webdev-biotool/flask/'+user_input



			#End cursor and connection to database
			c.close()
			conn.close()


			#Housekeeping
			gc.collect() #garbage collector for cleaning up unneeded stuff
			session['entered_id'] = True
			session['engaged'] = 'engaged' 
 

		return render_template('results.html', form=form, main_info = main_info, target_journals = target_journals)
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
