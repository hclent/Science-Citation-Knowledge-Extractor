from __future__ import print_function
from processors import *
import re, nltk, json, pickle, time
import json
from nltk.corpus import stopwords
import os.path
from multiprocessing import Pool
import logging

# source activate py34 #my conda python enviornment for this

#Create log
logging.basicConfig(filename='.multi_preprocess.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logging.info('Started')

#Stopwords
eng_stopwords = nltk.corpus.stopwords.words('english') #remove default english stopwords
logging.info('Stopword settings set')


#Set a PROCESSORS_SERVER environment variable.
#It may take a minute or so to load the large model files.
def connect_to_Processors(port_num):
  logging.warning('Connecting to the pyProcessors server may take a while')
  path = '/home/hclent/anaconda3/envs/pyProcessors/lib/python3.4/site-packages/py34/processors-server.jar'
  api = ProcessorsAPI(port=port_num, jar_path=path, keep_alive=True)
  logging.info('Connected to pyProcessors')
  return api


#Get all text docs for that pmid
#Returns list of strings e.g. ['1234_1.txt', '1234_2.txt', ...]
def retrieveDocs(pmid):
  docs = [] #list of strings
  folder = '/home/hclent/data/'
  files = os.listdir(folder)
  for f in files:
    if pmid in f and 'doc' not in f:
      docs.append(f) #str
  logging.debug('retrieved list of documents to processes')
  return docs



#Define function to call in parallel
#Annotation of docs will be the process in the pool
#Prints to JSON
def multiprocess(docs):
  pool = Pool(10)
  logging.debug('created 10 worker pools')
  t0 = time.time()
  results = pool.map_async(loadDocuments, docs) #docs = ['17347674_1.txt', '17347674_2.txt', '17347674_3.txt', ...]
  logging.debug('initialized map_async to loadDocs function with docs')
  #results = pool.map(loadDocuments, docs)
  #logging.debug('did map to loadDocs function with docs. NO async')     #no difference in performance between async and no async
  #pool.close()
  #logging.debug('closed pool')        #can exclude close and join and still have same error (not POSTing)
  #pool.join()
  #logging.debug('joined pool')
  print("pool work: done in %0.3fs." % (time.time() - t0))
  print(results.get())
  print(type(results))
  #print(results) #don't use .get() for map (without async)
  #no async returns 'list', with async returns 'multiprocessing.pool.MapResult'
  i = 0
  for biodoc in results.get():
    print(type(biodoc))
    save_path = '/home/hclent/data/'
    running_doc = docs[i]
    print(running_doc)
    pmid_i = running_doc.strip(".txt")
    completeName = os.path.join(save_path, ('doc_'+str(pmid_i)+'.json'))
    i += 1
    with open(completeName, 'w') as out:
      out.write(completeName.to_JSON())
      logging.debug('printed to json')
      print("* Dumped "+str(pmid_i) +" to JSON !!! ")
      print("\n")
  logging.info('done with everything!!')




#Input: String(text)
#Output: This method cleans the text of newline markups, DNA sequences, and some punctuation
def preProcessing(text):
  clean_text = re.sub('\\\\n', ' ', text) #replace \n with a space
  clean_text = re.sub('\([ATGC]*\)', '', clean_text) #delete any DNA seqs
  clean_text = re.sub('(\(|\)|\'|\]|\[|\\|\,)', '', clean_text) #delete certain punctuation
  clean_text = re.sub('\\\\xa0\d*\.?\d?[\,\-]?\d*\,?\d*', '', clean_text) #delete formatting around figures
  clean_text = re.sub('et\sal\.', ' ', clean_text) #replace "et al." with a space
  clean_text = re.sub('\s\d{4}[\s\.\,\;\-]?(\d{4})?', '', clean_text) #delete years
  clean_text = re.sub('[\.\,]\d{1,2}[\,\-]?(\d{1,2})?\,?', '', clean_text) #delete citations
  clean_text = re.sub('Fig\.|Figure', '', clean_text) #delete 'fig.' and 'figure'
  clean_text = clean_text.lower()
  logging.debug('cleaned the text')
  return clean_text



#Input: list of documents to process
#Output: Then it makes a "biodoc" using the PyProcessor's "BioNLP" Processor. This step takes a while for longer docs
#Output: This doc is saved to JSON.
def loadDocuments(doc):
  filenamePrefix = "/home/hclent/data/"
  filename = filenamePrefix + str(doc) #str(i)
  print("* Loading " +str(doc))
  logging.debug('found the the file '+str(doc))
  text = open(filename, 'r')
  text = text.read()
  logging.debug('read text with text.read()')
  clean_text = preProcessing(text)
  print("* cleaned " + str(doc) )
  print("* beginning annotation of "  + str(doc) )
  biodoc = api.bionlp.annotate(clean_text) #annotates to JSON
  print(type(biodoc))
  logging.debug('the biodoc of ' + str(doc) + ' is type ' + str(type(biodoc)))
  print("* " + str(doc) + " is type " + str(type(biodoc)))
  logging.debug('completed annotation of '  + str(doc) )
  return biodoc


api = connect_to_Processors(4242)
docs = retrieveDocs("17347674")
multiprocess(docs)
logging.info('Finished')


#################################################################

# def retrieveBioDocs(pmid):
#   bio_docs = [] #list of strings
#   folder = '/home/hclent/data/'
#   files = os.listdir(folder)
#   for f in files:
#     if pmid in f and 'doc' in f:
#       bio_docs.append(f) #str
#   logging.debug('retrieved list of documents to processes')
#   return bio_docs
#
#
# #Input: Processors annotated biodocs
# #Output: String of lemmas
# def grab_lemmas(biodoc):
#   lemmas_list = biodoc["lemmas"] #list
#   keep_lemmas = [w for w in lemmas_list if w.lower() not in eng_stopwords]
#   keep_lemmas = (' '.join(map(str, keep_lemmas))) #map to string. strings are necessary for the TFIDF
#   print(keep_lemmas)
#   return keep_lemmas
#
#
# #Input: Processors annotated biodocs
# #Output: List of named entities
# def grab_nes(biodoc):
#   ners_list = biodoc["nes"] #list
#   print(ners_list)
#   return ners_list
#
# #Input: Processors annotated biodocs (from JSON)
# #Output: List of strings of all lemmas
# def loadBioDoc(biodocs):
#   data_samples = []
#   nes_list = []
#   i = 1
#   filenamePrefix = '/home/hclent/data/'
#   for bd in biodocs:
#     filename = filenamePrefix + str(bd)
#     with open(filename) as data_file:
#       data = json.load(data_file)
#       print("loaded the json")
#       lemmas = grab_lemmas(data)
#       data_samples.append(lemmas)
#       nes = grab_nes(data)
#       nes_list.append(nes)
#     i +=1
#   return data_samples, nes_list


# bio_docs = retrieveBioDocs("17347674")
# print(bio_docs)
# ds, nl = loadBioDoc(bio_docs)
