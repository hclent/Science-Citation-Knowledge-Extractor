import re, operator, collections, json
from itertools import chain
from collections import defaultdict
from database_management import db_journals_and_dates, db_citations_retrieval, db_get_years_range


def flatten(listOfLists):
    return list(chain.from_iterable(listOfLists))

#gets journals and dates from the citations, and filters out repeat citaitons
def get_journals_and_dates(query, conn):
	pmcids = []
	journals = []
	dates = []
	pmid_list = query.split('+')  # list of string pmids
	for pmid in pmid_list:
		pmcid, js, ds = db_journals_and_dates(pmid, conn)
		pmcids.append(pmcid)
		journals.append(js)
		dates.append(ds)
	#flatten lists for usability
	all_pmcids = flatten(pmcids)
	all_journals = flatten(journals)
	all_dates = flatten(dates)
	# need to sift out duplicates!!!!
	combos = collections.OrderedDict.fromkeys(zip(all_pmcids, all_journals, all_dates))
	keep_journals = [ c[1] for c in combos]
	keep_dates = [c[2] for c in combos]
	return keep_journals, keep_dates



#this gets journal info for the bar chart in "statistics"
def statistics_dates_barchart(journals, dates, query, conn):
	x_vals = []
	y_vals = []

	yearDict = defaultdict(lambda:0)
	years_list = dates

	try:
		#first see if its in the db
		years_plus_range = db_get_years_range(query, conn) #use the years range of the entire query
		#print("got years range from db") #in db years_range looks like "2008+2017
		years_range = years_plus_range.split('+')
	except Exception as e: #if its not in the db, do it the manual way
		years_range = get_years_range(query)
		#print("calculate years_range by hand")


	# Associate journals with years
	journal_year = list(zip(journals, years_list)) #('Scientific Reports', '2016')

	for year in range(int(years_range[0]), int(years_range[1]) + 1):
		for pair in journal_year:
			if year == int(pair[1]):
				yearDict[year] +=1

	#print(yearDict)
	for year in sorted(yearDict):
		x_vals.append(year)
		y_vals.append(yearDict[year])

	return x_vals, y_vals


def journal_dates_barchart(journals, years_list, query, conn):
	x_vals = []
	y_vals = []

	yearDict = defaultdict(lambda: 0)

	years_range = get_years_range(query, conn)

	# Associate journals with years
	journal_year = list(zip(journals, years_list))  # ('Scientific Reports', '2016')

	for year in range(int(years_range[0]), int(years_range[1]) + 1):
		for pair in journal_year:
			if year == int(pair[1]):
				yearDict[year] += 1

	for year in sorted(yearDict):
		x_vals.append(year)
		y_vals.append(yearDict[year])

	return x_vals, y_vals


#Years range looks like (2008, 2017)
#The actual journals vis uses this once, but in the future its stored in the db and retrievable that way
def get_years_range(query, conn):
	years_list = []
	pmid_list = query.split('+') #list of string pmids
	#print(pmid_list)
	for pmid in pmid_list:
		apa_citations, db_journals, db_dates, db_urls = db_citations_retrieval(pmid, conn)

		#print(db_dates) #step1 : get years
		for d in db_dates:

			year = re.search('\d{4}', d)
			if year:
				y = int(year.group(0))
			else:
				try:
					y = int(d[0:4])
				except Exception as e:
					# we need an int no matter what...
					if d == "None":
						y = 2018 # The future... hopefully people will see thats impossible
					else:
						y = int(d[0:4])

			#if its ONLY numbers like 20170607, we're gonna gamble and just put the first 4.
			#Ugh some pe
			years_list.append(int(y)) #append year as an int

	sorted_years = (sorted(years_list))
	#This is just adjusting the years range so that the tickmarks on the journals vis look pretty
	if sorted_years[-1] - sorted_years[0] < 8: #if the years range is less than eight, need to adjust that
		start = sorted_years[0] - 8 #minus 9 from start year
		end = sorted_years[-1]
		years_range = (str(start), str(end))
	else:
		years_range = (str(sorted_years[0]), str(sorted_years[-1])) #define years range with (oldest, newest)
	return years_range


#Makes the json for the Journals visualization :)
def journals_vis(years_range, query, conn):
	journals, dates = get_journals_and_dates(query, conn)
	#print("JOURNALS VISUALIZATION")
	num_publications = len(journals) #UNIQUE publications only since duplicates have been flitered out
	#print("THERE ARE " + str(num_publications)+ " PUBLICATIONS")

	years_list = []
	#already changed to years by get_journals_and_dates above
	#print(dates) #step1 : get years
	#July 11, 2017: Updating the regex here to deal with DUMB PEOPLE WHO PUT 20170607 as the date in the metadata on pubmed
	for d in dates:

		year = re.search('\d{4}', d)
		if year:
			y = int(year.group(0))
		else:
			try:
				y = int(d[0:4])
			except Exception as e:
				if d == "None":
					y = 2018
				else:
					# we need an int no matter what...
					y = int(d[0:4])

		# y = re.sub('.[A-Z]{1}[a-z]{2}(.?\d{1,2})?', '', d) #delete month and day
		# y = re.sub('\D', '', y) #delete any extra letters that escaped for some reason
		years_list.append(y)
	#print(years_list)


	#Associate journals with years
	journal_year = list(zip(journals, years_list)) #('Scientific Reports', '2016')
	#print(journal_year)

	#Dictionary with "Journal": [year, year]
	#For looking up the years
	jyDict = defaultdict(list)
	i = 0
	for j in journals:
		if j == (journal_year[i][0]):
			jyDict[j] += [journal_year[i][1]]
			i+=1
	#len(jyDict) is the amount of UNIQUE journals
	number_unique_journals = len(jyDict)
	#print("IN "+str(number_unique_journals)+" UNIQUE JOURNALS")

	#Dictionary with "Journal": Number-of-publications
	#For looking up the total
	journalsTotalDict = defaultdict(lambda: 0)
	sum_j = 0
	for j in journals:
		journalsTotalDict[j] += 1
		sum_j +=1
	#print(journalsTotalDict)

	#Sorted by counts MOST --> LEAST
	sorted_dict = sorted(journalsTotalDict.items(), key=operator.itemgetter(1), reverse=True)
	unique_journals = [journal[0] for journal in sorted_dict]

	#NOT SORTED BY COUNTS
	#unique_journals = list(journalsTotalDict.keys())
	#print(unique_journals)

	publication_data = []
	for j in unique_journals:
		#print(j)
		#Initiate the dictionary for this journal
		journal_data = {
			"name": j,
			"articles": [], #[[year, number], [year, number]]
			"total": journalsTotalDict[j]   #total can get from journalsTotalDict with key (total is value)
		}
		#print("Years a paper was in this journal: "+ str(jyDict[j]))
		for year in range(int(years_range[0]), int(years_range[1]) + 1):
			#print("checking " +str(year) +" ...")
			sum_y = 0
			for entry in jyDict[j]:
				#print(" ... against "+str(entry))
				if year == int(entry):
					#print("The years match so I'm going to count now")
					sum_y+=1
				year_sum = [year, sum_y]
				#print(year_sum)
			journal_data["articles"].append(year_sum)

		publication_data.append(journal_data)

	range_info = [years_range, num_publications, number_unique_journals]
	#print(range_info)

	## add a TOTAL SUM row to the top of the journals visualization
	x, y = journal_dates_barchart(journals, years_list, query, conn)
	total_sum = sum(y)
	total_articles = [[year, count] for year, count in zip(x, y)]
	total_name = "TOTAL"
	top_row = {
		"name": total_name,
		"articles": total_articles,
		"total": total_sum
	}

	publication_data = [top_row] + publication_data

	#print("RANGE INFO: ")
	#Example range info: [('2008', '2016'), 165, 48] means years 2008-2016, 165 publications, 48 unique journals

	publication_data = json.dumps(publication_data)
	#print(publication_data)

	return publication_data, range_info





