from database_management import *

#total number
#total unique (overlap)


#pmcDict
#dict value [0] = num abstracts
#dict value [1] = num whole articles
#dict value [2] = num sentences
#dict value [3] = num tokens


def get_statistics(pmid_list):
    total = []
    unique_pmcids = []

    all_abstracts = []
    all_whole = []
    all_sents = []
    all_tokens = []

    for pmid in pmid_list:

        pmidDict, pmcDict = db_statistics(pmid)
        #print(pmidDict)
        total.append(pmidDict[pmid])
        print(pmcDict)
        for key, value in pmcDict.items():

            if key not in unique_pmcids:
                unique_pmcids.append(key)
            abstract = value[0]

            if abstract == 'yes':
                all_abstracts.append(abstract)

            whole = value[1]
            if abstract == 'yes':
                all_whole.append(whole)

            sent = value[2]
            all_sents.append(sent)

            token = value[3]
            all_tokens.append(token)


    sum_total = sum(total)
    unique = (len(unique_pmcids))
    sum_abstracts = len(all_abstracts)
    sum_whole = len(all_whole)
    sum_sents = sum(all_sents)
    sum_tokens = sum(all_tokens)
    statistics = [sum_total, unique, sum_abstracts, sum_whole, sum_sents, sum_tokens]
    print(statistics)
    return statistics




get_statistics(["18269575", "18952863"])





