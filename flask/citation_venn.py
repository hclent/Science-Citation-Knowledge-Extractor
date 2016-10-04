from database_management import db_citation_titles, db_citation_pmc_ids
from itertools import combinations
from collections import defaultdict


#pmid_list = ['18269575']
#pmid_list = ['18269575', '10467567', '9108111']

def make_venn(pmid_list):
    i = 0
    i_count = []

    all_pmcids = [] #[pmc1, pmc2, ... pmc-n]
    pmc_id_list = [] #list of lists [[pmid_ids paper1], [pmid_ids paper2]]

    citationsDict = defaultdict(list)

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
        pmid_Dict = {'sets': x, 'label': pmid_label, 'size': (len(pmc_ids))}
        #print(pmid_Dict)
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

            citationsDict[citation] = papers


        ### enumerate all possible dictionary values ####
        possible_vals = []
        s = ''
        combos = s.join(i_count) #format combos for 'combinations' method

        for n in range(2, len(pmid_list)):
            compare_all = (combinations(combos, int(n)))

            for c in compare_all:
                option_list = [] #dict keys are in list format [1, 2]
                #print(c)
                for number in c:
                    option_list.append(int(number))

                possible_vals.append(option_list)

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
        pass


    return venn_data

# venn_data = make_venn(pmid_list)
# print(venn_data)


############# GRAVEYARD #######################
# def format_venn(pmid_list):
#     sets = []
#     running_sets = []
#     all_titles = []
#     i_count = []
#     i = 0
#     for pmid in pmid_list:
#         titles = db_citation_titles(pmid)
#         x = [int(i)]
#         pmid_label = 'PMID:'+str(pmid)+' ('+str(len(titles))+')'
#         pmid_Dict = {'sets': x, 'label': pmid_label, 'size': (len(titles))}
#         print(pmid_Dict)
#         running_sets.append(pmid_Dict)
#         sets.append(pmid_Dict)
#         all_titles.append(titles)
#         i += 1
#
#     #if there is more than 1 pmid
#     if len(pmid_list) > 1:
#
#         n = len(running_sets)
#         for i in range(0, n):
#             i_count.append(str(i))
#         print(i_count)
#         s = ''
#         combos = s.join(i_count) #possible combination of elements
#         print(combos)
#
#         #pairwise
#         pairwise = (combinations(combos, 2))
#         for b in pairwise:
#             print(b)
#             both_titles = []
#             for t in all_titles[ int(b[0]) ]:
#                 if t in all_titles[ int(b[1]) ]:
#                     print(t)
#                     print("is in both papers")
#                     both_titles.append(t)
#             print(both_titles)
#             print(len(both_titles))
#             x = [int(b[0]), int(b[1])]
#             count_label = str(len(both_titles))
#             overlap_Dict = {'sets': x, 'label': count_label, 'size': int(len(both_titles))}
#             print(overlap_Dict)
#             sets.append(overlap_Dict)
#
#         #3-wise comparison
#         if n > 2:
#             threewise = (combinations(combos, 3))
#             for b in threewise:
#                 print(b)
#                 three_titles = []
#                 for t in all_titles[ int(b[0])  ]:
#                     if t in all_titles[int(b[1])] and t in all_titles[int(b[2])]:
#                         print(t)
#                         print("is in 3 papers")
#                         three_titles.append(t)
#                 print(three_titles)
#                 print(len(three_titles))
#                 x = [int(b[0]), int(b[1]), int(b[2])]
#                 count_label = str(len(three_titles))
#                 overlap_Dict = {'sets': x, 'label': count_label, 'size': int(len(three_titles))}
#                 print(overlap_Dict)
#                 sets.append(overlap_Dict)
#
#
#         #
#         #compare all
#         compareall = (combinations(combos, n))
#         for b in compareall:
#             print(b)
#             all_titles = []
#             for t in all_titles[ int(b[0]) ]:
#
#
#
#     print(sets)
#     for lilDicts in sets:
#         print(lilDicts)
#
#     #json needs double quotes, not single quotes
#     #print(sets)
#
#
# format_venn(pmid_list)


