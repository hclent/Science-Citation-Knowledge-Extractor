from __future__ import print_function
from processors import *
import re, nltk, json, pickle, time
import json
from nltk.corpus import stopwords
import os.path
from multiprocessing import Pool
import logging

# source activate py34 #my conda python environment for this

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
  api = ProcessorsAPI(port=port_num, jar_path=path, keep_alive=True, jvm_mem="-Xmx16G")
  logging.info('Connected to pyProcessors')
  rando_doc = api.bionlp.annotate("The mitochondria is the powerhouse of the cell.")
  logging.info('Annotated something random to initialize bioNLP Processor')
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
  pool = Pool(20)
  logging.debug('created worker pools')
  results = pool.map_async(loadDocuments, docs) #docs = ['17347674_1.txt', '17347674_2.txt', '17347674_3.txt', ...]
  logging.debug('initialized map_async to loadDocs function with docs')
  logging.debug('did map to loadDocs function with docs. WITH async')
  pool.close()
  logging.debug('closed pool')
  pool.join()
  logging.debug('joined pool')
  print(results.get())
  i = 0
  for biodoc in results.get():
    save_path = '/home/hclent/data/'
    running_doc = docs[i]
    pmid_i = running_doc.strip(".txt")
    completeName = os.path.join(save_path, ('doc_'+str(pmid_i)+'.json'))
    i += 1
    with open(completeName, 'w') as out:
      out.write(biodoc.to_JSON())
      logging.debug('printed to json')
  logging.info('Finished')




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
  logging.debug("NEW TASK")
  filenamePrefix = "/home/hclent/data/"
  filename = filenamePrefix + str(doc) #str(i)
  logging.debug('found the the file '+str(doc))
  text = open(filename, 'r')
  text = text.read()
  logging.debug('read text with text.read()')
  clean_text = preProcessing(text)
  logging.debug("* beginning annotation of "  + str(doc) )
  biodoc = api.bionlp.annotate(clean_text) #annotates to JSON  #thread safe?
  logging.debug('the biodoc of ' + str(doc) + ' is type ' + str(type(biodoc)))
  logging.debug("END OF TASK")
  return biodoc


#t1 = time.time()
api = connect_to_Processors(4343)
#docs = retrieveDocs("18269575") #"17347674"
#multiprocess(docs)
#print("Execute everything: done in %0.3fs." % (time.time() - t1))

###################################################
###################################
#Input: Processors annotated biodocs
#Output: String of lemmas
def retrieveBioDocs(pmid):
  biodocs = [] #list of strings
  folder = '/home/hclent/data/'
  files = os.listdir(folder)
  for f in files:
    if pmid in f and 'doc' in f:
      biodocs.append(f) #str
  logging.debug('retrieved list of Bio documents to work with')
  return biodocs


def grab_lemmas(biodoc):
  lemmas_list = biodoc.lemmas #list
  keep_lemmas = [w for w in lemmas_list if w.lower() not in eng_stopwords]
  keep_lemmas = (' '.join(map(str, keep_lemmas))) #map to string. strings are necessary for the TFIDF
  return keep_lemmas


#Input: Processors annotated biodocs
#Output: List of named entities
def grab_nes(biodoc):
  ners_list = biodoc.nes
  return ners_list

#Input: Processors annotated biodocs (from JSON)
#Output: List of strings of all lemmas
def loadBioDoc(biodocs):
  data_samples = []
  nes_list = []
  for bd in biodocs:
    filename = '/home/hclent/data/'+bd
    with open(filename) as jf:
      data = Document.load_from_JSON(json.load(jf))
      lemmas = grab_lemmas(data)
      data_samples.append(lemmas)
      nes = grab_nes(data)
      nes_list.append(nes)
  return data_samples, nes_list

#biodocs = retrieveBioDocs("18269575")
#data_samples, nes_list = loadBioDoc(biodocs)
