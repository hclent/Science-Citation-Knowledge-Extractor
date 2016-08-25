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
def do_kemeans(sparse_matrix):
    t0 = time.time()
    print("* Beginning k-means clustering ... ")
    num_clusters = 3
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
    i += 1

  #for 2D, just "go.Scatter"
  trace0 = go.Scatter3d(
    x = x0_coordinates,
    y = y0_coordinates,
    z = z0_coordinates,
    name = 'cluster 1',
    mode = 'markers',
    marker = dict(
      size = 10,
      color = 'rgba(152, 0, 0, .8)',
        line = dict(
            width = 2,
            color = 'rgb(0, 0, 0)'
        ),
        opacity=0.8
      )
    )

  trace1 = go.Scatter3d(
      x = x1_coordinates,
      y = y1_coordinates,
      z = z1_coordinates,
      name = 'cluster 2',
      mode = 'markers',
      marker = dict(
          size = 10,
          color = 'rgb(204, 204, 204)',
          line = dict(
              width = 2,
          ),
          opacity=0.8
      )
    )

  trace2 = go.Scatter3d(
    x = x2_coordinates,
    y = y2_coordinates,
    z = z2_coordinates,
    name = 'cluster 3',
    mode = 'markers',
    marker = dict(
        size = 10,
        color = 'rgba(156, 165, 196, 0.95)',
        line = dict(
            width = 2,
        ),
        opacity=0.8
      )
    )


  data = [trace0, trace1, trace2]
  layout = go.Layout(
    margin=dict(
        l=0,
        r=0,
        b=0,
        t=0
    )
  )

  #For 2D
  # layout = dict(title = 'K-means Clustering',
  #             yaxis = dict(zeroline = False),
  #             xaxis = dict(zeroline = False)
  #            )

  fig = dict(data=data, layout=layout)
  #plot_url = py.plot(data, filename='simple-3d-scatter') #print to online
  plot(data)
  print("done in %0.3fs." % (time.time() - t0))


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
clusters = do_kemeans(hX) #list of cluster assignments
coordinates = do_NMF(hX) #dimensionality reduction for visualization
plotKmeans(coordinates, clusters) #format for Plotly scatterplot





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