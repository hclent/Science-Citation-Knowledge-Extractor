import json, re, os
from fasttext import *
from operator import itemgetter

#Input: zipped results of fasttext.py
def embedding_json(results, query):
    json_out = {"nodes": [], "links": []}
    for r1 in results:
        g1 = int(r1[0])
        word1 = r1[1][0] #"p. vularis"
        formatted_word1 = str( re.sub("\s","_", word1) ) #needs to be formatted "p._vulgaris"
        json_out["nodes"].append({"id": formatted_word1, "group": (g1 + 1)}) #+1 so there's no Topic 0
        for r2 in results:
            g2 = int(r2[0])
            word2 = r2[1][0]
            formatted_word2 = str(re.sub("\s", "_", word2))
            if r1 == r2: #if they are exactly the same word, pass
                pass
            elif g1 == g2: #if the groups are the same, connect them.
                #print( word1 + ": " + str(g1) + " == " + word2 + ": " +str(g2)   )
                json_out["links"].append({"source": formatted_word1, "target": formatted_word2, "value": 1})
            else:
                pass #cosine sim data would be added here for words not in the same group

    #nodes actually need to be organized for the data vis
    ordered_nodes = sorted(json_out["nodes"], key=itemgetter('group'))
    links = json_out["links"]
    print_json = {"nodes": ordered_nodes, "links": links}
    save_path = '/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/fgraphs'  # in the folder of the last pmid
    completeName = os.path.join(save_path, ('fgraph_' + (str(query)) + '.json'))  # with the query for a name
    with open(completeName, "w") as outfile:
        json.dump(print_json, outfile)


# pmid_list = ['18952863','18269575']
# words, tags = get_words_tags(pmid_list) #list of words/tags per doc
# transformed_sentence = transform_text(words, tags)
# npDict = chooseTopNPs(transformed_sentence)
# top_nps= list(npDict.most_common(100))
# model = load_model('/home/hclent/tmp/fastText/16kmodel.vec')
# matrix = getNPvecs(top_nps, model)
# kmeans = KMeans(n_clusters=10, random_state=2).fit(matrix)
# res = list(zip(kmeans.labels_, top_nps))
# embedding_json(res)
# print("done!")

############ GRAVEYARD ##############
#reformat depreciated

#import pandas as pd

# data = pd.read_csv("/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_embed.csv", index_col=0)
# data.head()

# json_out = {"nodes": [], "links": []}
# for word, info in data.iterrows():
#     if word == "group":
#         continue
#     # print(word)
#     json_out["nodes"].append({"id": word, "group": int(info["group"])})
#     for target, weight in info.iteritems():
#         if target == "group":
#             continue
#
#         # print target, weight
#
#         if word == target:
#             # Case for sizing. TODO
#             pass
#         else:
#             # print target, weight
#             if weight == 1:
#                 # pass
#                 # print(target)
#                 # print(weight)
#                 json_out["links"].append({"source": word, "target": target, "value": 1})
#             if weight != 1:
#                 # print(target)
#                 # print(weight)
#                 json_out["links"].append({"source": word, "target": target, "value": weight})
#         #print("################")
#
#
# # print(json_out)
# with open("/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_embed.json", "w") as outfile:
#     json.dump(json_out, outfile)