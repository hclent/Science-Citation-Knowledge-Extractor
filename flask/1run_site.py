from flask import Flask, render_template, request, flash, url_for, redirect, session, g
from flask_wtf import Form
from wtforms import StringField, TextField, SelectField
import gc
import sqlite3
from content_management import runCrawler1
from content_management import runCrawler2
from content_management import runCosines
from cogeCrawled_db import connection

#If running in a virtualenv, must have modules also (pip) installed in that virtualenv! 
#flask, Flask-WTF, nltk, bs4, lxml, requests

app = Flask(__name__)

#Content from content_management.py
running_results = runCrawler1()
cosines = runCosines() 


#Main page: 
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
	return render_template('dashboard.html',running_results = running_results, cosines = cosines)



#Need to pass pmid in form to MainCrawler.py
class pmidForm(Form):
	pmid = TextField('PubmedID')


@app.route('/testing/', methods=["GET", "POST"])
def trying():
	form = pmidForm()
	if request.method == 'POST':
		flash('Posting: successful')
		pmid = form.pmid.data #pmid in test_page = pmid in pmidForm
		
		conn, c = connection()
		c.execute("INSERT INTO cogeCrawled (pmids) VALUES (?)", (pmid,)) #
		conn.commit()
		flash("Writing to database: successful")

		#End connection 
		c.close()
		conn.close()

		gc.collect() #garbage collector for cleaning up unneeded stuff
		session['entered_id'] = True
		session['pmid'] = pmid
	return render_template('dashboard.html', form=form)






#Handles 404 errors
@app.errorhandler(404)
def page_not_found(e):
	return("you shouldnt be here!!")



#Built in Flask debugger
if __name__ == '__main__':
	app.secret_key = 'super secret key'
	app.run(debug=True)

