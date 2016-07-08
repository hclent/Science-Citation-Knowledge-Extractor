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


#Create table
def create_table():
	c.execute('''CREATE TABLE IF NOT EXISTS cogeCrawled
		(post_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT, 
		datestamp TEXT, 
		pmids TEXT, 
		title TEXT,
		author TEXT,
		journal TEXT,
		comments TEXT)''') 


def data_entry():
	c.execute("INSERT INTO cogeCrawled VALUES(3, '05-25-2016', '43212', 'test2','author2', 'test2.com', 'PLOS')")
	conn.commit() #to save it to db
	
	c.execute("SELECT * FROM cogeCrawled")
	[print(row) for row in c.fetchall()]
	
	c.close()
	conn.close() 



def print_table():
	c.execute("SELECT * FROM cogeCrawled")
	[print(row) for row in c.fetchall()]


