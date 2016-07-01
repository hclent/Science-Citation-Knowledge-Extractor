from __future__ import print_function
from processors import *
import re, nltk, json, pickle, time
from json import JSONEncoder
from nltk.corpus import stopwords
import os.path

# source activate py34 #my conda python enviornment for this

# set a PROCESSORS_SERVER environment variable.
# It may take a minute or so to load the large model files.
#p = '/Users/hclent/anaconda3/envs/pyProcessors/lib/python3.4/site-packages/py34/processors-server.jar'
#api = ProcessorsAPI(port=8886, jar_path=p, keep_alive=True)


eng_stopwords = nltk.corpus.stopwords.words('english') #remove default english stopwords 
bio_stopwords = ['et', 'al', 'fig', 'author'] #add hand picked bio stopwords to stopwords
for word in bio_stopwords:
	eng_stopwords.append(word)


#Input: Data that you want to be JSONified
#Output: Data reformatted so it can be dumped to JSON
def dumper(obj):
  try:
    return obj.toJSON()
  except:
    return obj.__dict__

#Input: String(text), pmid, doc_num
#Output: This method cleans the text of newline markups, DNA sequences, and some punctuation
#Output: Then it makes a "biodoc" using the PyProcessor's "BioNLP" Processor. This step takes a while for longer docs
#Output: This doc is saved to JSON. 
#Output: pmid and doc_num are for naming the JSON filename 
def preProcessing(text, pmid, doc_num, api):
  print("* Preprocessing the text ... ")
  clean_text = re.sub('\\\\n', ' ', text) #replace \n with a space
  clean_text = re.sub('\([ATGC]*\)', '', clean_text) #delete any DNA seqs
  clean_text = re.sub('(\(|\)|\'|\]|\[|\\|\,)', '', clean_text) #delete certain punctuation
  clean_text = clean_text.lower()
  print("* Annotating with the Processors ...")
  print("* THIS WILL TAKE A WHILE ...")
  biodoc = api.bionlp.annotate(clean_text) #annotates to JSON
  print("* Successfully did the preprocessing !!!")
  print("* Dumping JSON ... ")
  save_path = '/Users/hclent/Desktop/webdev-biotool/flask/data/' #must save to static
  completeName = os.path.join(save_path, ('doc_'+(str(pmid))+'_'+str(doc_num)+'.json'))
  with open(completeName, 'w') as outfile:
    json.dump(biodoc, outfile, default=dumper, indent=2)
  print("* Dumped to JSON !!! ")


#Input: filehandle and max number of documents to process
#Output: JSONified annotated BioDoc 
def loadDocuments(maxNum, pmid, api):
  print("* Loading dataset...")
  i = 1
  filenamePrefix = "/Users/hclent/Desktop/webdev-biotool/flask/data/"+pmid+"_"
  print(filenamePrefix)
  for i in range(1, int(maxNum)+1):
    print("* Loading document #" + str(i) + " ...")
    filename = filenamePrefix + str(i) + ".txt"
    text = open(filename, 'r')
    text = text.read()
    preProcessing(text, pmid, i, api)
    i +=1
    print("\n")



###################################
#Input: Processors annotated biodocs
#Output: String of lemmas
def grab_lemmas(biodoc):
  lemmas_list = biodoc["lemmas"] #list 
  keep_lemmas = [w for w in lemmas_list if w.lower() not in eng_stopwords]
  keep_lemmas = (' '.join(map(str, keep_lemmas))) #map to string. strings are necessary for the TFIDF
  return keep_lemmas


#Input: Processors annotated biodocs
#Output: List of named entities 
def grab_nes(biodoc):
  ners_list = biodoc["nes"] #list 
  return ners_list

#Input: Processors annotated biodocs (from JSON)
#Output: List of strings of all lemmas 
def loadBioDoc(maxNum, pmid):
  data_samples = []
  nes_list = []
  i = 1
  filenamePrefix = '/Users/hclent/Desktop/webdev-biotool/flask/data/doc_'+(pmid)+'_'
  for i in range(1, maxNum+1):
    filename = filenamePrefix + str(i) + ".json"
    with open(filename) as data_file:
      data = json.load(data_file)
      lemmas = grab_lemmas(data)
      data_samples.append(lemmas)
      nes = grab_nes(data)
      nes_list.append(nes)
    i +=1
  return data_samples, nes_list

