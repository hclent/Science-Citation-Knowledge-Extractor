## test.py ###
import pandas as pd
import json

data = pd.read_csv("/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_embed.csv", index_col=0)
data.head()

json_out = {"nodes": [], "links": []}
for word, info in data.iterrows():
    if word == "group":
        continue
    # print(word)
    json_out["nodes"].append({"id": word, "group": int(info["group"])})
    for target, weight in info.iteritems():
        if target == "group":
            continue

        # print target, weight

        if word == target:
            # Case for sizing. TODO
            pass
        else:
            # print target, weight
            if weight == 1:
                # pass
                # print(target)
                # print(weight)
                json_out["links"].append({"source": word, "target": target, "value": 1})
            if weight != 1:
                # print(target)
                # print(weight)
                json_out["links"].append({"source": word, "target": target, "value": weight})
        #print("################")


# print(json_out)
with open("/home/hclent/repos/Webdev-for-bioNLP-lit-tool/flask/static/coge_embed.json", "w") as outfile:
    json.dump(json_out, outfile)