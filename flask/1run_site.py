from flask import Flask, render_template, request
from content_management import runCrawler1
from content_management import runCrawler2
from content_management import runCosines

#If running in the virtualenv, must have modules also installed in that virtualenv! 
#flask, nltk, bs4, lxml, requests

app = Flask(__name__)

running_results = runCrawler1() #Actually runs the program

crawl_results = runCrawler2() #Prints pickle obj

cosines = runCosines() #


@app.route('/cogecrawl')
def cogecrawl():
	return render_template('dashboard.html', cosines = cosines, crawl_results = crawl_results,
	running_results = running_results) #html file reference = actual dict


if __name__ == '__main__':
	app.run(debug=True)