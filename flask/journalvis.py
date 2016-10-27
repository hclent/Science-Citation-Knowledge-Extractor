import re
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

#get_years_range('18952863+18269575')


def journals_vis(journals, dates, years_range):
	print("JOURNALS VISUALIZATION CRAP")
	num_publications = len(journals)
	print("THERE ARE " + str(num_publications)+ " PUBLICATIONS")

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
	print("IN "+str(number_unique_journals)+" UNIQUE JOURNALS")

	#Dictionary with "Journal": Number-of-publications
	#For looking up the total
	journalsTotalDict = defaultdict(lambda: 0)
	sum = 0
	for j in journals:
		journalsTotalDict[j] += 1
		sum +=1
	#print(journalsTotalDict)
	unique_journals = list(journalsTotalDict.keys())
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
	print("RANGE INFO: ")
	print(range_info)
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


# d1 = ['2016 Jul 15', '2016 Jun 10', '2016 May 11', '2016 Jan 11', '2016 Mar 7', '2015 Nov 2', '2016 Apr 29', '2016 Apr 15', '2016 Mar 30', '2016 Mar 29', '2016 Mar 17', '2016 Feb 24', '2016 Feb 19', '2016 Jan 29', '2016 Jan 20', '2015 Dec 16', '2015 Nov 16', '2015 Oct 2', '2015 Aug 25', '2015 Jun 12', '2015 Jun 5', '2015 Apr 21', '2015 Mar 13', '2014 Oct 16', '2014 Nov 20', '2015 Feb 3', '2014 Dec 30', '2014 May 15', '2014 Nov 8', '2014 Oct 21', '2014 Oct 17', '2014 Oct 17', '2013 Oct 5', '2014 Jul 17', '2014 Jul 9', '2014 Jul 1', '2014 Jul 8', '2014 Jul 17', '2014 Aug 5', '2013 Dec 19', '2014 May 27', '2014 Apr 3', '2014 Mar 31', '2014 Mar 27', '2014 Feb 7', '2014 Mar 20', '2014 Mar 18', '2013 Dec 12', '2013 Dec 5', '2013 Oct 15', '2013 Oct 15', '2013 Sep 22', '2013 Nov 5', '2013 Aug 15', '2013 Jan 28', '2013 Jan 28', '2012 Nov 29', '2012 Aug 30', '2012 Oct 15', '2012 Oct 17', '2012 Sep 11', '2012 Sep 3', '2012 Aug 30', '2012 Aug 2', '2012 Jul 31', '2012 Jun 8', '2012 Jun 25', '2012 Jan 1', '2012 Apr 10', '2011 Mar 10', '2011 Jul 25', '2012 Jan 4', '2012 Feb 15', '2011 Dec 23', '2011 Dec 16', '2011 Nov 4', '2011 May 27', '2011 Sep 13', '2011 Apr 5', '2010 May 14', '2011 Mar 10', '2011 Feb 18', '2011 Feb 8', '2010 Jun 24', '2010 Jun 29', '2010 May 22', '2010 Feb 24', '2009 Oct 5', '2008 Dec', '2015 Aug 13']
# d2 = ['2016 Aug 28', '2016 May 6', '2016 Jul 25', '2016 Aug 17', '2016 Jun 10', '2016 Jun 10', '2016 May 16', '2016 Jun 29', '2016 Jan 11', '2016 May 17', '2016 Mar 29', '2016 Feb 1', '2016 Jan 29', '2015 Nov 11', '2015 Dec 16', '2015 Dec 16', '2015 Oct 2', '2015 Oct 7', '2015 Sep 11', '2015 Jun 5', '2015 Apr 9', '2014 Dec 2', '2014 Dec 30', '2014 Nov 27', '2014 Nov 8', '2014 Dec 11', '2014 Oct 17', '2013 Oct 5', '2014 Sep 11', '2014 Aug 25', '2014 Jun 30', '2014 Jul 9', '2014 Jul 8', '2014 Feb 20', '2014 Mar 28', '2013 May 21', '2013 Jun 3', '2013 Oct 15', '2013 Oct 15', '2013 Oct 29', '2013 Nov 11', '2013 Nov 5', '2013 Nov 1', '2013 Oct 31', '2013 Sep 13', '2013 Jul 2', '2013 Mar 14', '2013 Feb 5', '2013 Feb 1', '2013 Jan 4', '2012 Nov 23', '2012 Sep 3', '2012 Jul 31', '2012 Jun 8', '2012 Jun 25', '2012 May 24', '2012 Apr 10', '2012 Apr', '2011 Nov 17', '2011 Mar 11', '2010 May 24', '2011 Apr 18', '2011 Mar 24', '2011 Feb 25', '2010 Nov 23', '2010 Jun 29', '2010 May 1', '2010 May 22', '2010 May 15', '2008 Apr 24', '2009 Dec', '2009 Jul', '2009 Jun', '2008 Dec', '2008 Dec']
#
# j1 = ['G3: Genes|Genomes|Genetics', 'Plant Physiology', 'Chromosome Research', 'BMC Genomics', 'Molecular Genetics and Genomics', 'Nature genetics', 'PLoS ONE', 'Frontiers in Plant Science', 'Frontiers in Plant Science', 'Frontiers in Genetics', 'PLoS ONE', 'BioMed Research International', 'Journal of Experimental Botany', 'Frontiers in Plant Science', 'Scientific Reports', 'Frontiers in Plant Science', 'PLoS ONE', 'BMC Genomics', 'GigaScience', 'Frontiers in Plant Science', 'PLoS ONE', 'Frontiers in Plant Science', 'BMC Evolutionary Biology', 'Nucleic Acids Research', 'Nucleic Acids Research', 'Frontiers in Plant Science', 'BMC Plant Biology', 'Plant Molecular Biology Reporter / Ispmb', 'BMC Genomics', 'Plant Physiology', 'BMC Genomics', 'BMC Genomics', 'The Plant journal : for cell and molecular biology', 'The Plant Cell', 'BMC Genomics', 'BMC Genomics', 'GigaScience', 'PLoS Genetics', 'Philosophical Transactions of the Royal Society B: Biological Sciences', 'Molecular Biology and Evolution', 'Frontiers in Plant Science', 'International Journal of Molecular Sciences', 'Scientific Reports', 'PLoS ONE', 'Journal of Experimental Botany', 'PLoS ONE', 'PLoS ONE', 'International Journal of Molecular Sciences', 'BMC Genetics', 'BMC Bioinformatics', 'BMC Bioinformatics', 'BMC Evolutionary Biology', 'BMC Genomics', 'BMC Genomics', 'Genome Biology and Evolution', 'BMC Genomics', 'Nucleic Acids Research', 'Genome Biology', 'Journal of Experimental Botany', 'PLoS ONE', 'PLoS ONE', 'Bioinformatics', 'Frontiers in Plant Science', 'BMC Bioinformatics', 'Frontiers in Plant Science', 'The Plant Cell', 'BMC Bioinformatics', 'Mobile Genetic Elements', 'Proceedings of the National Academy of Sciences of the United States of America', 'Frontiers in Plant Science', 'Frontiers in plant science', 'Nucleic Acids Research', 'PLoS ONE', 'Plant Physiology', 'The Plant Cell', 'Genome Biology and Evolution', 'Genome Biology', 'The Plant Cell', 'Nucleic Acids Research', 'BMC Plant Biology', 'PLoS ONE', 'BMC Evolutionary Biology', 'BMC Plant Biology', 'Annals of Botany', 'PLoS Biology', 'Journal of Molecular Evolution', 'Bioinformatics', 'Genome Biology and Evolution', 'Genome Research', 'Plant Biotechnology Journal']
# j2 = ['Evolutionary Bioinformatics Online', 'Nucleic Acids Research', 'Proceedings of the National Academy of Sciences of the United States of America', 'Scientific Reports', 'Plant Physiology', 'Plant Physiology', 'Applied and Environmental Microbiology', 'BMC Evolutionary Biology', 'BMC Genomics', 'BMC Genomics', 'Frontiers in Genetics', 'Frontiers in Plant Science', 'Frontiers in Plant Science', 'Genome Biology and Evolution', 'Frontiers in Plant Science', 'Frontiers in Plant Science', 'BMC Genomics', 'Frontiers in Neuroscience', 'BMC Plant Biology', 'PLoS ONE', 'PeerJ', 'The Plant Cell', 'BMC Plant Biology', 'Plant and Cell Physiology', 'BMC Genomics', 'PLoS ONE', 'BMC Genomics', 'The Plant journal : for cell and molecular biology', 'G3: Genes|Genomes|Genetics', 'G3: Genes|Genomes|Genetics', 'BMC Genomics', 'BMC Genomics', 'GigaScience', 'The ISME Journal', 'The Plant Cell', 'Plant Biotechnology Reports', 'The New Phytologist', 'BMC Bioinformatics', 'BMC Bioinformatics', 'Genome Biology and Evolution', 'Proceedings of the National Academy of Sciences of the United States of America', 'BMC Genomics', 'PLoS ONE', 'PLoS Genetics', 'The Plant Cell', 'Frontiers in Plant Science', 'Genome Biology and Evolution', 'BMC Genomics', 'G3: Genes|Genomes|Genetics', 'Frontiers in Plant Science', 'Nucleic Acids Research', 'Bioinformatics', 'Frontiers in Plant Science', 'The Plant Cell', 'BMC Bioinformatics', 'PLoS ONE', 'Proceedings of the National Academy of Sciences of the United States of America', 'Genetics', 'Nucleic Acids Research', 'Nucleic Acids Research', 'BMC Plant Biology', 'BMC Bioinformatics', 'BMC Evolutionary Biology', 'PLoS ONE', 'Parasites & Vectors', 'PLoS Biology', 'BMC Evolutionary Biology', 'Journal of Molecular Evolution', 'Genes & Development', 'Nature', 'The Plant Cell', 'The Plant Cell', 'The Plant Cell', 'Plant Physiology', 'Genome Research']
#
# dates = d1+d2
# journals = j1+j2
# paper_dates_barchart(journals, dates, '18952863+18269575')

