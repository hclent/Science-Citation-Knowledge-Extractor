from flask import Flask, render_template, request, flash, url_for, redirect, session, g
from flask_wtf import Form
from wtforms import StringField, TextField, SelectField
import sqlite3
import gc
from content_management import runCrawler1
from content_management import runCrawler2
from content_management import runCosines
from cogeCrawled_db import connection
from MainCrawler import pmc_spider


#If running in a virtualenv, must have modules also (pip) installed in that virtualenv! 
#flask, Flask-WTF, nltk, bs4, lxml, requests

app = Flask(__name__)

#Content from content_management.py
pickle_results = runCrawler2()
cosines = runCosines() 


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

			if attempted_pmid == "1234":
				return redirect(url_for('cogecrawl'))

			else:
				error = "Invalid pmid. Try again."
	except Exception as e:
		#flash(e)
		return render_template("dashboard.html", error=error) #not showing error....
	return render_template('dashboard.html', pickle_results = pickle_results, cosines = cosines) #I should do the example page here :')



#Creat Form for handling user-entered pmid
#Need to pass pmid in form to MainCrawler.py
class pmidForm(Form):
	pmid = TextField('PubmedID')



#Getting Flask-WTFs to work with sqlite3 here
#This function uses user entry to run the MainCrawler.py 
#Prints the results under the "Information" tab
#User entered pmid is entered into sqlite3 database
#Problem: Can't do lists with sqlite... 
@app.route('/results/', methods=["GET", "POST"])
def trying():
	form = pmidForm()
	if request.method == 'POST':
		flash('Posting to form: successful')

		
		pmid = form.pmid.data #THIS IS THE USER INPUT FROM THE FORM
	

		#Using userinput to start crawl
		user_input = str(pmid)
		titles, urls, authors = pmc_spider(1, user_input) #returns 3 lists #only crawling first page!
		main_info = list(zip(titles, authors, urls))


		#sqlite3 database entries
		conn, c = connection()
		c.execute("INSERT INTO cogeCrawled (pmids) VALUES (?)", (pmid,)) #put user pmid into db
		conn.commit()
		flash("Writing PubmedID to database: successful")


		#End connection 
		c.close()
		conn.close()


		#Housekeeping
		gc.collect() #garbage collector for cleaning up unneeded stuff
		session['entered_id'] = True
		session['pmid'] = pmid

	return render_template('results.html', form=form, main_info = main_info)




#Handles 404 errors
@app.errorhandler(404)
def page_not_found(e):
	return("you shouldnt be here!!")



#Built in Flask debugger
if __name__ == '__main__':
	app.secret_key = 'super secret key'
	app.run(debug=True)

