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
  t0 = time.time()
  results = pool.map_async(loadDocuments, docs) #docs = ['17347674_1.txt', '17347674_2.txt', '17347674_3.txt', ...]
  logging.debug('initialized map_async to loadDocs function with docs')
  print("pool established: done in %0.3fs." % (time.time() - t0))
  logging.debug('did map to loadDocs function with docs. WITH async')
  pool.close()
  logging.debug('closed pool')
  pool.join()
  logging.debug('joined pool')
  print(results.get()) #timeout=1, returns None for all
  i = 0
  for biodoc in results.get():
    save_path = '/home/hclent/data/'
    running_doc = docs[i]
    print("Preparing to dump "+ str(running_doc) + " to JSON")
    pmid_i = running_doc.strip(".txt")
    completeName = os.path.join(save_path, ('doc_'+str(pmid_i)+'.json'))
    i += 1
    with open(completeName, 'w') as out:
      out.write(biodoc.to_JSON())
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
  logging.debug("NEW TASK")
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
  logging.debug("* beginning annotation of "  + str(doc) )
  biodoc = api.bionlp.annotate(clean_text) #annotates to JSON  #thread safe?
  logging.debug('the biodoc of ' + str(doc) + ' is type ' + str(type(biodoc)))
  print("* " + str(doc)   + " is type " + str(type(biodoc)))
  logging.debug("END OF TASK")
  return biodoc


t1 = time.time()
api = connect_to_Processors(4343)
docs = retrieveDocs("18269575") #"17347674"
multiprocess(docs)
logging.info('Finished')
print("Execute everything: done in %0.3fs." % (time.time() - t1))

