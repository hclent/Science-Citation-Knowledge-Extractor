from flask import Flask, render_template, request, flash, url_for, redirect
from content_management import runCrawler1
from content_management import runCrawler2
from content_management import runCosines

#If running in the virtualenv, must have modules also installed in that virtualenv! 
#flask, nltk, bs4, lxml, requests

app = Flask(__name__)

#Content:
running_results = runCrawler1() #Actually runs the program
cosines = runCosines() #


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



@app.errorhandler(404)
def page_not_found(e):
	return("you shouldnt be here!!")


#Built in Flask debugger
if __name__ == '__main__':
	app.secret_key = 'super secret key'
	app.run(debug=True)