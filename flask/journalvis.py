import re, operator
from collections import defaultdict
from database_management import db_citations_retrieval


def get_years_range(query):
	years_list = []
	pmid_list = query.split('+') #list of string pmids
	#print(pmid_list)
	for pmid in pmid_list:
		apa_citations, db_journals, db_dates, db_urls = db_citations_retrieval(pmid)
		#print(db_dates) #step1 : get years
		for d in db_dates:
			y = re.sub('.[A-Z]{1}[a-z]{2}(.?\d{1,2})?', '', d) #delete month and day
			y = re.sub('\D', '', y) #delete any extra letters that escaped for some reason
			years_list.append(int(y)) #append year as an int
		#print("------------------------")
	sorted_years = (sorted(years_list))
	if sorted_years[-1] - sorted_years[0] < 8: #if the years range is less than eight, need to adjust that
		start = sorted_years[0] - 8 #minus 9 from start year
		end = sorted_years[-1]
		years_range = (str(start), str(end))
	else:
		years_range = (str(sorted_years[0]), str(sorted_years[-1])) #define years range with (oldest, newest)
	#print("YEARS RANGE: " + str(years_range))

	return years_range


def journals_vis(journals, dates, years_range):
	#print("JOURNALS VISUALIZATION CRAP")
	num_publications = len(journals)
	#print("THERE ARE " + str(num_publications)+ " PUBLICATIONS")

	years_list = []

	#print(dates) #step1 : get years
	for d in dates:
		y = re.sub('.[A-Z]{1}[a-z]{2}(.?\d{1,2})?', '', d) #delete month and day
		y = re.sub('\D', '', y) #delete any extra letters that escaped for some reason
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
	sum = 0
	for j in journals:
		journalsTotalDict[j] += 1
		sum +=1
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
			sum = 0
			for entry in jyDict[j]:
				#print(" ... against "+str(entry))
				if year == int(entry):
					#print("The years match so I'm going to count now")
					sum+=1
				year_sum = [year, sum]
				#print(year_sum)
			journal_data["articles"].append(year_sum)

		publication_data.append(journal_data)

	range_info = [years_range, num_publications, number_unique_journals]
	#print("RANGE INFO: ")
	#print(range_info)
	#Get some info about the publication before changing it to a string for the json
	#Year range, number of publications, number of unique journals
	publication_data = re.sub('\'', '\"', str(publication_data)) #json needs double quotes, not single quotes
	return (publication_data, range_info)



def paper_dates_barchart(journals, dates, query):
	x_vals = []
	y_vals = []


	yearDict = defaultdict(lambda:0)
	years_list = []

	years_range = get_years_range(query)

	for d in dates:
		y = re.sub('.[A-Z]{1}[a-z]{2}(.?\d{1,2})?', '', d) #delete month and day
		y = re.sub('\D', '', y) #delete any extra letters that escaped for some reason
		years_list.append(y)

	# Associate journals with years
	journal_year = list(zip(journals, years_list)) #('Scientific Reports', '2016')

	for year in range(int(years_range[0]), int(years_range[1]) + 1):
		for pair in journal_year:
			if year == int(pair[1]):
				yearDict[year] +=1

	for year in sorted(yearDict):
		x_vals.append(year)
		y_vals.append(yearDict[year])

	return x_vals, y_vals

