from database_management import db_citation_titles
from itertools import combinations
import re


#pmid_list = ['18269575', '18952863']
pmid_list = ['18269575', '18952863', '10467567', '9108111', '15504555']



def format_venn(pmid_list):
    sets = []
    running_sets = []
    all_titles = []
    i_count = []
    i = 0
    for pmid in pmid_list:
        titles = db_citation_titles(pmid)
        x = [int(i)]
        pmid_label = 'PMID:'+str(pmid)+' ('+str(len(titles))+')'
        pmid_Dict = {'sets': x, 'label': pmid_label, 'size': (len(titles))}
        print(pmid_Dict)
        running_sets.append(pmid_Dict)
        sets.append(pmid_Dict)
        all_titles.append(titles)
        i += 1

    #if there is more than 1 pmid
    if len(pmid_list) > 1:

        n = len(running_sets)
        for i in range(0, n):
            i_count.append(str(i))
        print(i_count)
        s = ''
        combos = s.join(i_count) #possible combination of elements
        print(combos)

        #pairwise
        pairwise = (combinations(combos, 2))
        for b in pairwise:
            print(b)
            both_titles = []
            for t in all_titles[ int(b[0]) ]:
                if t in all_titles[ int(b[1]) ]:
                    print(t)
                    print("is in both papers")
                    both_titles.append(t)
            print(both_titles)
            print(len(both_titles))
            x = [int(b[0]), int(b[1])]
            count_label = str(len(both_titles))
            overlap_Dict = {'sets': x, 'label': count_label, 'size': int(len(both_titles))}
            print(overlap_Dict)
            sets.append(overlap_Dict)

        #3-wise comparison
        if n > 2:
            threewise = (combinations(combos, 3))
            for b in threewise:
                print(b)
                three_titles = []
                for t in all_titles[ int(b[0])  ]:
                    if t in all_titles[int(b[1])] and t in all_titles[int(b[2])]:
                        print(t)
                        print("is in 3 papers")
                        three_titles.append(t)
                print(three_titles)
                print(len(three_titles))
                x = [int(b[0]), int(b[1]), int(b[2])]
                count_label = str(len(three_titles))
                overlap_Dict = {'sets': x, 'label': count_label, 'size': int(len(three_titles))}
                print(overlap_Dict)
                sets.append(overlap_Dict)


        #
        #compare all
        compareall = (combinations(combos, n))
        for b in compareall:
            print(b)
            all_titles = []
            for t in all_titles[ int(b[0]) ]:
                


    print(sets)
    for lilDicts in sets:
        print(lilDicts)

    #json needs double quotes, not single quotes
    #print(sets)




format_venn(pmid_list)
