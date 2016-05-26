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
#missing urls column because i forgot a column like a fool!!!!! 


# def data_entry():
# 	c.execute("INSERT INTO cogeCrawled VALUES(3, '05-25-2016', '43212', 'test2','author2', 'test2.com', 'PLOS')")
# 	conn.commit() #to save it to db
	
# 	c.execute("SELECT * FROM cogeCrawled")
# 	[print(row) for row in c.fetchall()]
	
# 	c.close()
# 	conn.close() 


# def dynamic_data_entry():
# 	date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
# 	c.execute("INSERT INTO cogeCrawled (datestamp, keyword, value) VALUES(?, ?, ?, ?)", #sqlite uses ? instead of %s
# 		(unix, date, keyword, value)) 
	#conn.commit()


def print_table():
	c.execute("SELECT * FROM cogeCrawled")
	[print(row) for row in c.fetchall()]


#connection()
#create_table()

#print_table()



##################################################################################      
# pmids_info.db as of May 23-2016
# entries 4, and on from website :') 

# post_id     datestamp   pmids       title       author       journal      comments  
# ----------  ----------  ----------  ----------  -----------  -----------  ----------
# 1           05-23-2016  1234        test title  test author  testurl.com  nature    
# 2           05-25-2016  4321        test2       author2      test2.com    PLOS      
# 3           05-25-2016  43212       test2       author2      test2.com    PLOS      
# 4                       zzzz                                                        
# 5                       wwwww                                                       
# 6                       alksjrlkje                                                  
# 7                       abcdefg                                                     
# 8                       2249436                                                     
# 9                       1643513                                                     
# 10                      8353730                                                     
# 11                      17045299                                                    
# 12                      19194716  