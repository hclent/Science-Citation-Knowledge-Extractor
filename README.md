# webdev-for-bioNLP-tool
* WORK IN PROGRESS

BASICS
--------------------------------------
* This will be the web interface for the bioNLP lit tool
* Built with Python3, Flask, and templates from getboostrap
* In order to run this site you must first activate a virtualenv (". venv/source/activate" or "conda create") in the main directory
* Once virtual enviornment is running in your terminal, pip install flask, Flask-WTF, nltk, bs4, lxml, requests
* Then from the flask directory run "python3 1run_site.py"
* Open browser and go to "http://localhost:5000/cogecrawl" to see the site and some example data
* Enter a PubmedID into the textbook and hit "Submit", then click "Information" to see the corresponding data for that PubmedID

ORGANIZATION
--------------------------------------
* The file content_management.py imports functions from other code from the bioNLP tool to return the information which I want on the site.
* Currently, only some of the bioNLP literature tool is in the "Flask" folder. More to come.
* Pickle objects will not be part of the final code, but are currently there as placeholders to save time. 

* 1run_site.py currently prints results from Maincrawler.py for authors, pub name, etc. to the site
* Also prints cosine similarity scores of two genetics texts (science1.txt and science2.txt) compared to a movie review (movie.txt).

* my_sqlitedb.py is in progress. It is for managing pmids entered on the site and the corresponding results from Maincrawler.py

TEMPLATES
--------------------------------------
See dashboard.html for main page content and some example results (dashboard.html is a continuation of header.html)
See results.html for the results of the user queried pubmedID. Currently only printing the titles, authors, and urls of citations
