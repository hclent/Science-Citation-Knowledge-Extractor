from cosineClassify import loadMessages
from cosineClassify import cosineSimilarityScore


############# COSINE STUFF to print in dashboard.html ##############
def runCosines():
	doc1 = "/Users/hclent/Desktop/webdev-biotool/flask/science1.txt"
	doc2 = "/Users/hclent/Desktop/webdev-biotool/flask/science2.txt"
	movie_review = "/Users/hclent/Desktop/webdev-biotool/flask/movie.txt"

	vecs1 = loadMessages(doc1)
	vecs2 = loadMessages(doc2)
	vecs3 = loadMessages(movie_review)

	cosines = cosineSimilarityScore(vecs1, vecs2, vecs3)

	return cosines

