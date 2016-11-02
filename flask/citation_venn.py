from database_management import db_citation_titles, db_citation_pmc_ids
from itertools import combinations
from collections import defaultdict


#pmid_list = ['18269575']
#pmid_list = ['18269575', '18952863']

def make_venn(pmid_list):
    i = 0
    i_count = []

    all_pmcids = [] #[pmc1, pmc2, ... pmc-n]
    pmc_id_list = [] #list of lists [[pmid_ids paper1], [pmid_ids paper2]]

    citationsDict = defaultdict(list)
    possible_vals = []

    venn_data = []
    for pmid in pmid_list:
        #add pmc_ids to master list for making citationsDict
        pmc_ids = db_citation_pmc_ids(pmid)
        pmc_id_list.append(pmc_ids)
        for id in pmc_ids:
            all_pmcids.append(id)

        #do first entries for venn
        x = [int(i)]
        pmid_label = 'PMID:'+str(pmid)+' ('+str(len(pmc_ids))+')'
        #href_base = 'https://www.ncbi.nlm.nih.gov/pubmed/'+str(pmid)
        #href_label = str('<a href="'+href_base+'">'+str(pmid_label)+'</a>')
        pmid_Dict = {'sets': x, 'label': pmid_label, 'size': (len(pmc_ids))}
        #print(pmid_Dict)
        venn_data.append(pmid_Dict)
        i_count.append(str(i))

        i+=1



    #if there is more than one pmid,
    if (len(pmid_list)) > 1:
        #print("more than 1")
        #make dictionary of {citating: ['papers', 'that', 'are', 'cited'] }
        for citation in all_pmcids:
            j = 0
            papers = []

            for pmclist in pmc_id_list:
                if citation in pmclist:
                    #print(str(citation) + " is in P"+str(j))
                    if j not in papers:
                        papers.append(j)
                j +=1
                #print(papers)
            citationsDict[citation] = papers


        ### enumerate all possible dictionary values ####
        s = ''
        combos = s.join(i_count) #format combos for 'combinations' method
        #print(combos)
        for n in range(2, len(pmid_list)+1):
            compare_all = (combinations(combos, int(n)))

            for c in compare_all:
                option_list = [] #dict keys are in list format [1, 2]
                #print("c")
                #print(c)
                for number in c:
                    #print("number")
                    #print(number)
                    option_list.append(int(number))

                possible_vals.append(option_list)

        #print("inside possible_vals")
        #print(possible_vals)

    #print(citationsDict)

    #print("OUTSIDE possible vales")
    #print(possible_vals)

    #if there were multiple pmids, must look at overlap
    if citationsDict:
        for p_vals in possible_vals:
            #print(p_vals)
            sum = 0
            for key, value in citationsDict.items():
                if value == p_vals:
                    #print(str(key)+' has '+str(value)+' == '+str(p_vals))
                    sum += 1
            #print(sum)
            if sum > 0:
                count_label = str(sum)
                overlap_Dict = {'sets': p_vals, 'label': count_label, 'size': sum }
                #print(overlap_Dict)
                venn_data.append(overlap_Dict)


    #if there was only 1 input pmid, no other work needed
    if not citationsDict:
        #print("there's no citationDict")
        pass


    return venn_data

#venn_data = make_venn(pmid_list)
#print(venn_data)



