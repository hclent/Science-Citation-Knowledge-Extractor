import json, re, os, logging
from operator import itemgetter

logging.basicConfig(filename='.app.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

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
    logging.info(completeName)
    with open(completeName, "w") as outfile:
        json.dump(print_json, outfile)


