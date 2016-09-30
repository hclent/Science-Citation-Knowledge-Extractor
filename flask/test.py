from multi_preprocess import *



def blahblah(biodocs):
  total_sentences = []
  total_tokens = []
  for bd in biodocs: #e.g. doc_26502977_3.json
    pmid_i = bd.strip('.json')
    pmid = re.sub('(doc\_)', '', pmid_i)
    pmid = re.sub('\_\d*', '', pmid)
    filename = '/home/hclent/data/'+(str(pmid))+'/'+bd
    with open(filename) as jf:
      data = Document.load_from_JSON(json.load(jf))
      num_sentences = data.size
      print(num_sentences)
      total_sentences.append(num_sentences)
      for i in range(0, num_sentences):
        s = data.sentences[i]
        num_tokens = s.length
        print(num_tokens)
        total_tokens.append(num_tokens)

  sum = 0
  for sent in total_sentences:
    sum += sent
  print("TOTAL SENTS: " + str(sum))

  total = 0
  for tokens in total_tokens:
    total += tokens
  print("TOTAL TOKENS: " + str(total))





biodocs = retrieveBioDocs("18952863")
blahblah(biodocs)



