from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.decomposition import NMF
from sklearn.cluster import KMeans
import sys, pickle, math, random, numpy, time
import plotly.plotly as py
import plotly.graph_objs as go
py.sign_in('hclent', 'eeg49e9880')
from plotly.offline import plot
from multi_preprocess import * #mine


#Input: Eata_samples (list of lists containing strings)
#Output: Sparse matrix, l2 normalization for preserving Euclidean distance
def get_hashing(data):
  t0 = time.time()
  print("* Making hashing vectorizor with the data ...")
  hasher = HashingVectorizer(stop_words='english', ngram_range=(1,3), norm='l2', non_negative=True) #l2 projected on the euclidean unit sphere
  hX = hasher.fit_transform(data)
  print("done in %0.3fs." % (time.time() - t0))
  return hX, hasher


#Input: High dimensional (sparse) matrix
#Output: Clusters
# labels = km.labels_
# centroids = km.cluster_centers_
def do_kemeans(sparse_matrix, k_clusters):
    t0 = time.time()
    print("* Beginning k-means clustering ... ")
    num_clusters = int(k_clusters)
    km = KMeans(init='k-means++', n_clusters=num_clusters)
    km.fit(sparse_matrix)
    clusters = km.labels_.tolist()
    print("done in %0.3fs." % (time.time() - t0))
    return clusters


#Non-Negative Matrix Factorization
#Input: sparse matrix
#Output: list of Cartesian coordinates for each document vector
def do_NMF(sparse_matrix):
  t0 = time.time()
  print("* Performing NMF on sparse matrix ... ")
  nmf = NMF(n_components=3)
  coordinates = nmf.fit_transform(sparse_matrix)
  print("done in %0.3fs." % (time.time() - t0))
  return(coordinates)


#Function for making a 3D plot in Plotly
#Input: Cartesian coordinates and document cluster assignments
#Output: 3D scatter plot
def plotKmeans(coordinates, clusters):
  t0 = time.time()
  print("* Preparing to plot now ... ")
  x0_coordinates = []
  y0_coordinates = []
  z0_coordinates = []
  x1_coordinates = []
  y1_coordinates = []
  z1_coordinates = []
  x2_coordinates = []
  y2_coordinates = []
  z2_coordinates = []
  x3_coordinates = []
  y3_coordinates = []
  z3_coordinates = []
  x4_coordinates = []
  y4_coordinates = []
  z4_coordinates = []
  i = 0
  for vectors in coordinates:
    if clusters[i] == 0:
      x0_coordinates.append(vectors[0])
      y0_coordinates.append(vectors[1])
      z0_coordinates.append(vectors[2])
    if clusters[i] == 1:
      x1_coordinates.append(vectors[0])
      y1_coordinates.append(vectors[1])
      z1_coordinates.append(vectors[2])
    if clusters[i] == 2:
      x2_coordinates.append(vectors[0])
      y2_coordinates.append(vectors[1])
      z2_coordinates.append(vectors[2])
    if clusters[i] == 3:
      x3_coordinates.append(vectors[0])
      y3_coordinates.append(vectors[1])
      z3_coordinates.append(vectors[2])
    if clusters[i] == 4:
      x4_coordinates.append(vectors[0])
      y4_coordinates.append(vectors[1])
      z4_coordinates.append(vectors[2])
    i += 1


  print("done in %0.3fs." % (time.time() - t0))
  return(x0_coordinates, y0_coordinates, z0_coordinates,
         x1_coordinates, y1_coordinates, z1_coordinates,
         z2_coordinates, z2_coordinates, z2_coordinates,
         z3_coordinates, z3_coordinates, z3_coordinates,
         z4_coordinates, z4_coordinates, z4_coordinates)


bdocs1 = retrieveBioDocs("18269575")
data_samples, neslist1 = loadBioDoc(bdocs1)
bdocs2 = retrieveBioDocs("18952863")
datas2, neslist2 = loadBioDoc(bdocs2)

for docs in datas2:
    data_samples.append(docs)
print("added datas 2 to data_samples")

hX, hasher = get_hashing(data_samples)
print(hX.toarray())
print(hX.shape)
print()
clusters = do_kemeans(hX, 5) #list of cluster assignments
coordinates = do_NMF(hX) #dimensionality reduction for visualization
x0_coordinates, y0_coordinates, z0_coordinates, x1_coordinates, y1_coordinates, z1_coordinates, z2_coordinates, z2_coordinates, z2_coordinates,z3_coordinates, z3_coordinates, z3_coordinates,z4_coordinates, z4_coordinates, z4_coordinates = plotKmeans(coordinates, clusters) #format for Plotly scatterplot


print(z4_coordinates)


####### GRAVEYARD ##########
#Hashing Vector is better than TF-IDF for K-means because of Eucliean distance

# from sklearn.feature_extraction.text import TfidfVectorizer
# from sklearn.feature_extraction.text import TfidfTransformer
# from sklearn.decomposition import TruncatedSVD
# from sklearn.pipeline import make_pipeline
# from sklearn.preprocessing import Normalizer

# def get_tfidf(data): #data should be a list of strings for the documents
#   print("* Making tfidf with the data ...")
#   tfidf_vectorizer = TfidfVectorizer(stop_words='english', ngram_range=(1, 3), norm='l2') #l2 projected on the euclidean unit sphere
#   tfidf = tfidf_vectorizer.fit_transform(data)
#   return tfidf, tfidf_vectorizer

# tfidfX, tfidf_vectorizer = get_tfidf(data_samples)
# print(tfidfX.toarray())
# print(tfidfX.shape)
# print()

# #Truncated SVD (LSA) for dimensionality reduction
# #For plotting only (don't want to give Dimen-Reduced data to kmeans!)
# def dimen_reduce(sparse_matrix):
#   print("* Performing SVD on sparse matrix ... ")
#   svd = TruncatedSVD(n_components=3, n_iter=100)
#   normalizer = Normalizer(copy=False)
#   lsa = make_pipeline(svd, normalizer)
#   X = lsa.fit_transform(sparse_matrix)
#   explained_variance = svd.explained_variance_ratio_.sum()
#   print("Explained variance of the SVD step: {}%".format(
#     int(explained_variance * 100)))
#   return X

# svdX = dimen_reduce(hX)
# print(svdX)
# print(svdX.shape) #(84, 3)