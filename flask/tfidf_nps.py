from multi_preprocess import retrieveBioDocs
from processors import *
from fasttext import flatten
from sklearn.feature_extraction.text import CountVectorizer, TfidfVectorizer
from collections import Counter

# Results seem to be distributed alphabetically? w,x,y,z words have higher TF-IDF score
# a,b,c,d words have lower TF-IDF score
# For extracting the valuable NPs with TF-IDF, n-grams give MUCH better results

def doc_words_tags(pmid_list):
    words = []
    tags = []
    for pmid in pmid_list:
        biodocs = retrieveBioDocs(pmid)
        #print(len(biodocs))
        for doc in biodocs:
            doc_words = []
            doc_tags = []
            with open(doc) as jf:
                data = Document.load_from_JSON(json.load(jf))
                num_sentences = data.size
                for i in range(0, num_sentences):
                    s = data.sentences[i]
                    s_words = s.lemmas
                    s_tags = s.tags
                    doc_words.append(s_words)
                    doc_tags.append(s_tags)
            doc_words = flatten(doc_words)
            doc_tags = flatten(doc_tags)
            words.append(doc_words)
            tags.append(doc_tags)
    return words, tags

def transform_text(words, tags):
    transformed_tokens = []
    np_stack = []
    join_char = " " #instead of ("_")
    for i, w in enumerate(words):
        if tags[i].startswith("NN"):
            np_stack.append(w)
        else:
            if len(np_stack) > 0:
                np = join_char.join(w.lower() for w in np_stack)
                transformed_tokens.append(np)
                # reset stack
                np_stack = []
            # append current token
            transformed_tokens.append(w)
    return transformed_tokens


# pmid_list = ['18952863', '18269575']
# #pmid_list = ['9108111', '10467567',  '23226128', '22073781']
# words, tags = doc_words_tags(pmid_list)
# corpus = []
# for i in range(0, len(words)):
#     flat_words = words[i]
#     flat_tags = tags[i]
#     transformed = transform_text(flat_words,flat_tags) #list
#     doc = [t for t in  transformed if len(t.split(' ')) > 1 ]
#     doc2string = str((' ').join(doc))
#     corpus.append(doc2string)
# vec = TfidfVectorizer(min_df=5, ngram_range=(2,3))
# counts = vec.fit_transform(corpus)
# counts.todense()
# tfidf = vec.vocabulary_
# top = Counter(tfidf)#.most_common(1000)
# print(top)


# these need to be clustered somehow, not with kmeans
