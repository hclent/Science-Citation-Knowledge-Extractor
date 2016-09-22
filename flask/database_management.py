import sqlite3
import time
import datetime
import random

#Basic SQLITE3 structure
#Database in webdev-biotool for managing pmids and scraped webpages
conn = sqlite3.connect(database='pmids_info.db') #connect to database
c = conn.cursor() #cursor


#Define connection and cursor
def connection():
	conn = sqlite3.connect(database='pmids_info.db') #connect to database
	c = conn.cursor() #cursor
	return conn, c


#Create table for inputPapers
def create_table_input():
	c.execute('''CREATE TABLE IF NOT EXISTS inputPapers
		(post_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
		datestamp TEXT,
		pmid TEXT,
		title TEXT,
		author TEXT,
		journal TEXT,
		pubdate TEXT,
		url TEXT)''')


#Create table for citations of the inputPapers
def create_table_citations():
	c.execute('''CREATE TABLE IF NOT EXISTS citations
		(post_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
		datestamp TEXT,
		pmcid TEXT,
		title TEXT,
		author TEXT,
		journal TEXT,
		pubdate TEXT,
		citesPmid TEXT,
		url TEXT)''')



def test_data_entry():
	c.execute("INSERT INTO citations VALUES(0, '09-21-2016', '000', 'title','author', 'journal', 'pubdate', 'pmid', 'www.website.come')")
	conn.commit() #to save it to db
	
	c.execute("SELECT * FROM citations")
	[print(row) for row in c.fetchall()]
	
	c.close()
	conn.close()


def print_inputPapers():
	c.execute("SELECT * FROM inputPapers")
	[print(row) for row in c.fetchall()]


def print_citations():
	c.execute("SELECT * FROM citations")
	[print(row) for row in c.fetchall()]


########### TO BE PHASED OUT ###################
#the cogeCrawled table will become obsolete
# def create_table():
# 	c.execute('''CREATE TABLE IF NOT EXISTS cogeCrawled
# 		(post_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
# 		datestamp TEXT,
# 		pmids TEXT,
# 		title TEXT,
# 		author TEXT,
# 		journal TEXT,
# 		comments TEXT)''')


# def print_cogeCrawled():
# 	c.execute("SELECT * FROM cogeCrawled")
# 	[print(row) for row in c.fetchall()]