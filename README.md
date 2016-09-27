# webdev-for-bioNLP-tool
* WORK IN PROGRESS

BASICS
--------------------------------------
* This will be the web interface for the bioNLP lit tool
* Built with Python3, Flask, and templates from getbootstrap
* See requirements.txt
* Displays default information about CoGe publications.
* Enter a PubmedID into the textbox and hit "Submit", then explore the tabs to see the corresponding data for that PubmedID


ORGANIZATION
--------------------------------------
* The file content_management.py imports functions from other code from the bioNLP tool to return the information which I want on the site.
* Currently, only some of the BioNLP Literature Tool is in the "Flask" folder. More to come.
* "Journals" tab visualizes citing publications in journals per year
* "NEs" tab visualises Named Entities from the citing papers, with a word cloud and heatmap
* "Topics" tab visualizes two topic models
* "Clustering" tab visualizes k-means clustering
* "Citations" lists citations with hyperlinks to the papers
* "SciFi" is in progress. 


TEMPLATES
--------------------------------------
See dashboard.html for main page content and some example results (dashboard.html is a continuation of header.html)
See results.html for the results of the user queried pubmedID.
