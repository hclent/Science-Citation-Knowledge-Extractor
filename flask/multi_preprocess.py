from __future__ import print_function
from processors import *
import re, nltk, json, pickle, time
import os.path
from multiprocessing import Pool
import logging
from configapp import app
from database_management import db_citation_pmc_ids, pmcidAnnotated  #mine

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
  path = (app.config['PATH_TO_JAR'])
  api = ProcessorsAPI(port=port_num, jar_path=path, keep_alive=True, jvm_mem="-Xmx10G")
  logging.info('Connected to pyProcessors')
  rando_doc = api.bionlp.annotate("The mitochondria is the powerhouse of the cell.")
  logging.info('Annotated something random to initialize bioNLP Processor')
  return api


#Initialize API as global variable in multi_preprocess.py
#We will keep it as a global variable because initializing the bioNLP processor has a cost
#Do not want to initialize the processors over and over
api = connect_to_Processors(4343)
#for odin path, lets do port 4342


#Get all text files for that pmid
#Returns list of dictionaries [{"pmcid": 1234, "filename": /path/to/file}, {}, {}, ...]
#Only retrieve txts for documents that have not been annotated before

def retrieveDocs(pmid, conn):
  docs = [] #list of strings
  #connect to db for citing id's
  db_pmcids = db_citation_pmc_ids(pmid, conn)
  for pmcid in db_pmcids:
    annotated = pmcidAnnotated(pmcid, conn)
    if annotated == 'yes':
      pass
    if annotated == 'no':
      pass
    if annotated == 'empty':
      prefix = pmcid[0:3]
      suffix = pmcid[3:6]

      #~path/123/456/12345678.txt
      folder = str(prefix)+'/'+str(suffix)+'/'+str(pmcid)+'.txt' #
      filename = os.path.join((app.config['PATH_TO_CACHE']), folder)
      docdict = {"pmcid": pmcid, "filepath": filename}
      docs.append(docdict)
  logging.debug('retrieved list of documents to processes')
  return docs



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
  #a doc is a dict {"pmcid": 1234, "filename": /path/to/file/name}
  pmcid = doc["pmcid"]
  filepath = doc["filepath"]

  #logging.debug("NEW TASK")
  #api = connect_to_Processors(4343) #could connect each time if don't want a global var
  logging.debug('found the the file '+str(filepath))
  #read the text
  text = open(filepath, 'r')
  text = text.read()
  logging.debug('read text with text.read()')
  clean_text = preProcessing(text)
  logging.debug("* beginning annotation of "  + str(pmcid) )
  biodoc = api.bionlp.annotate(clean_text) #annotates to JSON  #thread safe?
  logging.debug('the biodoc of ' + str(pmcid) + ' is type ' + str(type(biodoc)))

  ### let's try some trouble shooting for documents that failed
  if "processors.ds.Document" not in str(type(biodoc)):
    logging.debug("annotating this document failed! :( Let's try again HALF SIZE .... ")
    logging.debug(clean_text)
    doc_length = len(clean_text)
    half_length = int(doc_length * 0.5)
    half_clean_text = clean_text[0:half_length]
    biodoc = api.bionlp.annotate(half_clean_text)
    logging.debug('SECOND TRY: the biodoc of ' + str(pmcid) + ' is type ' + str(type(biodoc)))
    #if it fails AGAIN
    if "processors.ds.Document" not in str(type(biodoc)):
      #cut into a quarter size
      doc_length = len(clean_text)
      qt_length = int(doc_length * 0.25)
      qt_clean_text = clean_text[0:qt_length]
      biodoc =  api.bionlp.annotate(qt_clean_text)
      logging.debug('THIRD TRY (QUARTER SIZE): the biodoc of ' + str(pmcid) + ' is type ' + str(type(biodoc)))
      #if that still fails, return a blank json dict
      if "processors.ds.Document" not in str(type(biodoc)):
        fake_clean_text = 'error annotating document'
        biodoc = api.bionlp.annotate(fake_clean_text)

  logging.debug("END OF TASK")
  return biodoc, pmcid

#Define function to call in parallel
#Annotation of docs will be the process in the pool
#Prints to JSON
def multiprocess(docs):
  if len(docs) == 0: #if there are no docs, abort
    logging.info("no documents to annotate")
    pass
  else:
    t1 = time.time()
    pool = Pool(75)
    logging.debug('created worker pools')
    results = pool.map_async(loadDocuments, docs) #docs = ['17347674_1.txt', '17347674_2.txt', '17347674_3.txt', ...]
    logging.debug('initialized map_async to loadDocs function with docs')
    logging.debug('did map to loadDocs function with docs. WITH async')
    pool.close()
    logging.debug('closed pool')
    pool.join()
    logging.debug('joined pool')
    logging.info(results.get())

    for biodoc, pmcid in results.get():
      prefix = pmcid[0:3]
      suffix = pmcid[3:6]
      save_path = str(prefix)+'/'+str(suffix)+'/'+str(pmcid)+'.json' #look in folder that matches pmcid
      completeName = os.path.join((app.config['PATH_TO_CACHE']), save_path)

      with open(completeName, 'w') as out:
        out.write(biodoc.to_JSON())
        logging.debug('printed to json')

    logging.info("All biodoc creations: done in %0.3fs." % (time.time() - t1))
    logging.info('Finished')

###################################################
#BIODOC HANDLING

# Input: Pmid
# Output: list of dictionaries for all annotated citing pmcids ["pmcid": 123, ]
def retrieveBioDocs(pmid, conn):
  #print("retrieving biodocs")
  biodocs = [] #list of strings

  db_pmcids = db_citation_pmc_ids(pmid, conn)
  for pmcid in db_pmcids:
    annotated = pmcidAnnotated(pmcid, conn)
    if annotated == 'no':
      pass
    if annotated == 'empty':
      pass
    if annotated == 'yes':
      prefix = pmcid[0:3]
      suffix = pmcid[3:6]

      folder = str(prefix) + '/' + str(suffix) + '/' + str(pmcid) + '.json'  #
      filename = os.path.join((app.config['PATH_TO_CACHE']), folder)

      biodict = {"pmcid": pmcid, "jsonpath": filename}
      biodocs.append(biodict)
  logging.debug('retrieved list of Bio documents to work with')
  return biodocs


def grab_lemmas_and_tags(biodoc):
  lemmas_list = biodoc.lemmas #list
  tags_list = biodoc.tags
  lemmas_with_tags = list(zip(lemmas_list, tags_list))
  keep = [lt for lt in lemmas_with_tags if lt[0].lower() not in eng_stopwords]

  keep_lemmas = [lt[0] for lt in keep]
  keep_tags = [lt[1] for lt in keep]
  return keep_lemmas, keep_tags


#Input: Processors annotated biodocs
#Output: List of named entities
def grab_nes(biodoc):
  ners_list = biodoc.nes
  return ners_list



#Input: Processors annotated biodocs (from JSON)
#Output: list of dicts containing {pmcid, lemmas, nes, num_sentences, num_tokens}
def loadBioDoc(biodocs):
  t1 = time.time()

  loadedBioDocs = []

  for doc in biodocs:
    pmcid = doc["pmcid"]
    jsonpath = doc["jsonpath"]

    #IMPORTANT NOTE: MAY 10, 2017: "data_samples" being replaced with "lemmas" for clarity!!!!
    #biodict = {"pmcid": pmcid, "lemmas": [], "nes": [], "num_sentences": [], "num_tokens": [], "tags": []}
    biodict = {"pmcid": pmcid, "jsonpath": jsonpath, "nes": [], "lemmas": []} # "tags": [] gets added

    token_count_list = []

    with open(jsonpath) as jf:
      data = Document.load_from_JSON(json.load(jf))
      #print(type(data)) is <class 'processors.ds.Document'>
      num_sentences = data.size
      #biodict["num_sentences"].append(num_sentences)
      biodict["num_sentences"] = num_sentences
      for i in range(0, num_sentences):
        s = data.sentences[i]
        num_tokens = s.length
        token_count_list.append(num_tokens)

      num_tokens = sum(token_count_list)
      #biodict["num_tokens"].append(num_tokens)
      biodict["num_tokens"] = num_tokens

      lemmas, tags = grab_lemmas_and_tags(data)
      biodict["lemmas"].append(lemmas) #lemmas is a LIST
      # biodict["tags"].append(tags) #tags is a LIST
      biodict["tags"] = tags

      nes = grab_nes(data)
      biodict["nes"].append(nes)

      loadedBioDocs.append(biodict)

  logging.info("Done assembling lemmas, nes, sent+token counts: done in %0.3fs." % (time.time() - t1))
  return loadedBioDocs

