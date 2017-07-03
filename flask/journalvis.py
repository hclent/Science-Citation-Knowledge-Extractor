import re, operator, collections, json
from itertools import chain
from collections import defaultdict
from database_management import db_citations_retrieval, db_journals_and_dates


def flatten(listOfLists):
    return list(chain.from_iterable(listOfLists))

#gets journals and dates from the citations, and filters out repeat citaitons
def get_journals_and_dates(query):
	pmcids = []
	journals = []
	dates = []
	pmid_list = query.split('+')  # list of string pmids
	for pmid in pmid_list:
		pmcid, js, ds = db_journals_and_dates(pmid)
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
#TODO: make stacked barchart for each paper in query
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

#Years range looks like (2008, 2017)
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


def journals_vis(years_range, query):
	journals, dates = get_journals_and_dates(query)


	#print("JOURNALS VISUALIZATION")
	num_publications = len(journals) #UNIQUE publications only since duplicates have been flitered out
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

	## add a TOTAL SUM row to the top of the journals visualization
	x, y = paper_dates_barchart(journals, dates, query)
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
	#Get some info about the publication before changing it to a string for the json
	#Year range, number of publications, number of unique journals
	#publication_data = re.sub('\'', '\"', str(publication_data)) #json needs double quotes, not single quotes
	publication_data = json.dump(publication_data)

	return (publication_data, range_info)


#
# years_range = get_years_range("18952863+18269575")
# publication_data, range_info = journals_vis(years_range, "18952863+18269575") #[('2008', '2016'), 165, 48]
# print(publication_data)
# print(range_info)
#


