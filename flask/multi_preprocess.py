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
  path = '/home/hclent/anaconda3/envs/py34/lib/python3.4/site-packages/processors/processors-server.jar'
  api = ProcessorsAPI(port=port_num, jar_path=path, keep_alive=True, jvm_mem="-Xmx25G")
  logging.info('Connected to pyProcessors')
  rando_doc = api.bionlp.annotate("The mitochondria is the powerhouse of the cell.")
  logging.info('Annotated something random to initialize bioNLP Processor')
  return api


#Initialize API as global variable in multi_preprocess.py
#We will keep it as a global variable because initializing the bioNLP processor has a cost
#Do not want to initialize the processors over and over
api = connect_to_Processors(4343)


#Get all text docs for that pmid
#Returns list of strings e.g. ['1234_1.txt', '1234_2.txt', ...]
def retrieveDocs(pmid):
  docs = [] #list of strings
  folder = '/home/hclent/data/'+pmid+'/' #look in folder named after pmid
  files = os.listdir(folder)
  for f in files:
    if pmid in f and 'doc' not in f:
      #dont read .json or pickle objs
      if '.txt' in f and 'self' not in f: #don't include 'self' docs
        docs.append(f) #str
  logging.debug('retrieved list of documents to processes')
  return docs


#Define function to call in parallel
#Annotation of docs will be the process in the pool
#Prints to JSON
def multiprocess(docs):
  t1 = time.time()
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
    running_doc = docs[i]
    pmid_i = running_doc.strip(".txt")
    pmid = re.sub('\_\d*', '', pmid_i)
    save_path = '/home/hclent/data/'+str(pmid)+'/' #save in folder named after pmid
    completeName = os.path.join(save_path, ('doc_'+str(pmid_i)+'.json'))
    i += 1
    with open(completeName, 'w') as out:
      out.write(biodoc.to_JSON())
      logging.debug('printed to json')
  logging.info("All biodoc creations: done in %0.3fs." % (time.time() - t1))
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
  #api = connect_to_Processors(4343) #could connect each time if don't want a global var
  pmid_i = doc.strip(".txt")
  pmid = re.sub('\_\d*', '', pmid_i)
  filenamePrefix = "/home/hclent/data/"+str(pmid)+'/'
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

###################################################

#Input: Processors annotated biodocs
#Output: String of lemmas
##### THESE ARE NOT SORTED AND SHOULD PROBABLY BE SORTED ####
def retrieveBioDocs(pmid):
  print("retrieving biodocs")
  biodocs = [] #list of strings
  folder = '/home/hclent/data/'+pmid+'/'
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
  t1 = time.time()
  data_samples = []
  nes_list = []

  total_sentences = []

  total_tokens = []
  sum_tokens = []


  bd = biodocs[0] #just grab an example to get file name info

  pmid_i = bd.strip('.json')
  pmid = re.sub('(doc\_)', '', pmid_i)
  pmid = re.sub('\_\d*', '', pmid)


  #need to load the files in ORDER so that data_samples is in order
  for b in range(1, len(biodocs)+1):
    #filename = '/home/hclent/data/'+(str(pmid))+'/'+bd
    filename = '/home/hclent/data/' + (str(pmid)) + '/' + "doc_"+str(pmid)+"_"+str(b)+".json"
    #print(filename)
    doc_tokens= []
    with open(filename) as jf:
      data = Document.load_from_JSON(json.load(jf))
      #print(type(data)) is <class 'processors.ds.Document'>
      num_sentences = data.size
      total_sentences.append(num_sentences)
      for i in range(0, num_sentences):
        s = data.sentences[i]
        num_tokens = s.length
        doc_tokens.append(num_tokens)
      lemmas = grab_lemmas(data)
      data_samples.append(lemmas)
      nes = grab_nes(data)
      nes_list.append(nes)
      total_tokens.append(doc_tokens)
  logging.info("Done assembling lemmas and nes: done in %0.3fs." % (time.time() - t1))

  #add up tokens
  for sents in total_tokens:
    sum = 0
    for tokens in sents:
      sum += tokens
    sum_tokens.append(sum)

  logging.info("Done assembling sent counts and token counts")
  return data_samples, nes_list, total_sentences, sum_tokens

# docs = retrieveDocs("18269575")
# print(docs)
# multiprocess(docs)
# biodocs = retrieveBioDocs("18269575")
# loadBioDoc(biodocs)
# data_samples, nes_list, total_sentences, sum_tokens = loadBioDoc(biodocs)
# print(total_sentences)
# print(sum_tokens)