import sqlite3
import time
import datetime
import random
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import style
style.use('fivethirtyeight')


#Basic SQLITE3 structure
#To be implemented into webdev-biotool for managing pmids and scraped webpages


#Define connection and Cursor
conn = sqlite3.connect('mydatabase.db') #database
c = conn.cursor() #cursor


#Create a table
def create_table():
	c.execute('CREATE TABLE IF NOT EXISTS stuffToPlot(unix REAL, datestamp TEXT, keyword TEXT, value REAL)') #all caps for pure sql #normal for non-sql


#Enter data into table by hand
def data_entry():
	c.execute("INSERT INTO stuffToPlot VALUES(14523425234, '05-17-2016', 'Python', 8)")
	conn.commit() #to save it to db
	c.close()
	conn.close() 


#Automatically populate table with some random data
def dynamic_data_entry():
	unix = time.time()
	date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H: %M: %S'))
	keyword = 'Python'
	value = random.randrange(0,10)
	c.execute("INSERT INTO stuffToPlot (unix, datestamp, keyword, value) VALUES(?, ?, ?, ?)", #sqlite uses ? instead of %s
		(unix, date, keyword, value)) 
	conn.commit()


#Reading from database
def read_from_db():
	c.execute("SELECT * FROM stuffToPlot") #'Select *' selects everything with the cursor
	#data = c.fetchall() 
	#print(data)
	for row in c.fetchall(): 
		print(row[0])


#Graphing information from data
def graph_data():
	c.execute('SELECT unix, value FROM stuffToPlot')
	dates = []
	values = []
	for row in c.fetchall():
		dates.append(datetime.datetime.fromtimestamp(row[0]))
		values.append(row[1])
	plt.plot_date(dates, values, '-')
	plt.show()




# create_table()
# data_entry()


# for i in range(10):
# 	dynamic_data_entry()
# 	time.sleep(1)


#read_from_db()


#graph_data()



c.close()
conn.close()