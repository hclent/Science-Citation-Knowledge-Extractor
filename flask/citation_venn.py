from database_management import db_citation_pmc_ids
from itertools import combinations
from collections import defaultdict


def make_venn(pmid_list, conn):
    i = 0
    i_count = []

    all_pmcids = [] #[pmc1, pmc2, ... pmc-n]
    pmc_id_list = [] #list of lists [[pmid_ids paper1], [pmid_ids paper2]]

    citationsDict = defaultdict(list)
    possible_vals = []

    venn_data = []
    for pmid in pmid_list:
        #add pmc_ids to master list for making citationsDict
        pmc_ids = db_citation_pmc_ids(pmid, conn)
        pmc_id_list.append(pmc_ids)
        for id in pmc_ids:
            all_pmcids.append(id)

        #do first entries for venn
        x = [int(i)]
        pmid_label = 'PMID:'+str(pmid)+' ('+str(len(pmc_ids))+')'

        pmid_Dict = {'sets': x, 'label': pmid_label, 'size': (len(pmc_ids))}

        venn_data.append(pmid_Dict)
        i_count.append(str(i))

        i+=1

    #if there is more than one pmid,
    if (len(pmid_list)) > 1:
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
                for number in c:
                    option_list.append(int(number))

                possible_vals.append(option_list)


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
                venn_data.append(overlap_Dict)


    #if there was only 1 input pmid, no other work needed
    if not citationsDict:
        #print("there's no citationDict")
        pass


    return venn_data
