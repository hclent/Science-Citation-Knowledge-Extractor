import pickle
from content_management import db_inputPapers_retrieval


def statsSelfInfo(query):
    input_click_citations = []
    pmid_list = query.split('+')  # list of string pmids
    for user_input in pmid_list:
        apa = db_inputPapers_retrieval(user_input)
        url = "https://www.ncbi.nlm.nih.gov/pubmed/"+str(user_input)
        href_label = str('<a href="' + url + '">' + str(apa) + '</a>')
        input_click_citations.append(href_label)
    return(input_click_citations)

query = "18952863+18269575"
input_click_citations = statsSelfInfo(query)
print(input_click_citations)


