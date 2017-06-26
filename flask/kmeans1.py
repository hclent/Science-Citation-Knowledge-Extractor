from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.decomposition import NMF
from sklearn.cluster import KMeans
import pickle, random, time, logging


logging.basicConfig(filename='.app.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

#Input: Eata_samples (list of lists containing strings)
#Output: Sparse matrix, l2 normalization for preserving Euclidean distance
def get_hashing(data):
  list_data_strings = [' '.join(map(str, d)) for d in data] #map to string. strings are necessary for the TFIDF
  t0 = time.time()
  logging.info("* Making hashing vectorizor with the data ...")
  hasher = HashingVectorizer(stop_words='english', ngram_range=(1,3), norm='l2', non_negative=True)
  #l2 projected on the euclidean unit sphere
  hX = hasher.fit_transform(list_data_strings)
  logging.info("done in %0.3fs." % (time.time() - t0))
  return hX, hasher


#Input: High dimensional (sparse) matrix
#Output: Clusters
# labels = km.labels_
# centroids = km.cluster_centers_
# n_jobs for speeding up processing.... go back and figure out how to get it to work :/
def do_kemeans(sparse_matrix, k_clusters):
    t0 = time.time()
    logging.info("* Beginning k-means clustering ... ")
    num_clusters = int(k_clusters)
    km = KMeans(init='k-means++', n_clusters=num_clusters, n_jobs=10)
    km.fit(sparse_matrix)
    clusters = km.labels_.tolist()
    logging.info("done in %0.3fs." % (time.time() - t0))
    return clusters


#Non-Negative Matrix Factorization
#Input: sparse matrix
#Output: list of Cartesian coordinates for each document vector
def do_NMF(sparse_matrix):
  t0 = time.time()
  logging.info("* Performing NMF on sparse matrix ... ")
  nmf = NMF(n_components=3)
  coordinates = nmf.fit_transform(sparse_matrix)
  logging.info("done in %0.3fs." % (time.time() - t0))
  return(coordinates)


#Function for making a 3D plot in Plotly
#Input: Cartesian coordinates and document cluster assignments
#Output: 3D scatter plot
#coordinates must be a zip where ([vector], 'title')
def plotKmeans(coordinates, clusters):
  t0 = time.time()
  logging.info("* Preparing to plot now ... ")
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

  titles0 = []
  titles1 = []
  titles2 = []
  titles3 = []
  titles4 = []

  i = 0
  for vectors in coordinates:
    if clusters[i] == 0:
      x0_coordinates.append(vectors[0][0])
      y0_coordinates.append(vectors[0][1])
      z0_coordinates.append(vectors[0][2])
      titles0.append(vectors[1])
    if clusters[i] == 1:
      x1_coordinates.append(vectors[0][0])
      y1_coordinates.append(vectors[0][1])
      z1_coordinates.append(vectors[0][2])
      titles1.append(vectors[1])
    if clusters[i] == 2:
      x2_coordinates.append(vectors[0][0])
      y2_coordinates.append(vectors[0][1])
      z2_coordinates.append(vectors[0][2])
      titles2.append(vectors[1])
    if clusters[i] == 3:
      x3_coordinates.append(vectors[0][0])
      y3_coordinates.append(vectors[0][1])
      z3_coordinates.append(vectors[0][2])
      titles3.append(vectors[1])
    if clusters[i] == 4:
      x4_coordinates.append(vectors[0][0])
      y4_coordinates.append(vectors[0][1])
      z4_coordinates.append(vectors[0][2])
      titles4.append(vectors[1])
    i += 1


  logging.info("done in %0.3fs." % (time.time() - t0))
  return(x0_coordinates, y0_coordinates, z0_coordinates,
         x1_coordinates, y1_coordinates, z1_coordinates,
         x2_coordinates, y2_coordinates, z2_coordinates,
         x3_coordinates, y3_coordinates, z3_coordinates,
         x4_coordinates, y4_coordinates, z4_coordinates,
         titles0, titles1, titles2, titles3, titles4)

