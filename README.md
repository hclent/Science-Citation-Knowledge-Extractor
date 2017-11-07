# Science Citation Knowledge Extractor (SCKE)


## About SCKE
* SCKE is an open source tool that helps biomedical researchers understand how their work is being used by others, by analyzing the content in papers that cite them. This tool uses natural language processing and machine learning to extract the prominent themes and concepts discussed in the citing documents. By seeing what kinds of topics are discussed by citing articles, researchers can better understand how their work influences their peers and various disciplines in science. Additionally, SCKE allows biomedical researchers to explore other statistics about the publications that cite them, such as where citations are published (Journals), the distribution of keywords (Keywords), the similarity of papers to each other (Clustering), the similarity of papers to other famous works (TextCompare), and general statistics about the citations (Statistics).
* Built with `Python3` and `Flask` using `Biopython`, `py-processors`, `NLTK`, `Scikit-Learn`, `Numpy`, `Gensim`, `Fasttext`, `Plotly`, `D3`, and more.
* SCKE proudly leverages the BioNLP Processor created by the Computational Language Understanding Lab. See [here](https://github.com/clulab/processors) for original code. See [here](https://github.com/myedibleenso/py-processors) for python wrappper library, `py-processors`.


## Instructions for using the site
* Video to come :)

## Running SCKE
#### Docker container for running the site:
* Docker container(s) to come!

#### Instructions for setting up the site without a Docker container

See our wiki for [installation instructions](https://github.com/hclent/Science-Citation-Knowledge-Extractor/wiki/Installation-Instructions)!

## Features

### Journals
Where are your citers publishing their work? Could you be influencing authors with publications in Nature?

![Top of Journals vis](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/static/images/screenies/journals.png)

### Keywords
Visualizations the key words in your papers with either a wordcloud, heatmap, or clustermap! Sort by categories of keywords: Bioprocess, Cell-lines, Cellular components, Family, Gene or gene products, Organs, Simple chemicals, Sites, Species, Tissue-types.

#### wordcloud

![Wordcloud for Species and Families](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/static/images/screenies/wordcloud.png)

#### heatmap

![Heatmap for Cell Components](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/static/images/screenies/heatmap.png)

### Topics
See the latent themes in the documents that site you! Check out `Embeddings`, which uses word vectors to cluster important words into topics. If you have a large number of citations, `LSA` (Latent Semantic Analysis) and `LDA` (Latend Dirichlet Analysis) might work well for you too.

![Topics generated with Embeddings](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/static/images/screenies/embedding.png)

### Clustering
Project citing documents into 3D space! The most similar documents will be close together in the graph or the same color. Zoom and click to explore.

![Document Clustering](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/static/images/screenies/kmeans.png)

### Citations
A list of your citations

![Top of List of citations](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/static/images/screenies/citation.png)

### Statistics
Get the big picture about your citers -- How many (unique) papers cite you? How many of them are unique? How many does SCKE take into account? Which years do you have the most citers for your paper(s)?

![Statistics](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/static/images/screenies/statistics.png)

### TextCompare
See how similar your citing publications are to famous works like "On the Origin of Species" or "Sherlock Holmes". If your publication(s) are available, see how similar your citers are to your work!

![TextCompare results against input document](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/static/images/screenies/textcompare.png)

**Note:** The yellow bar in this visualization is one of the input papers (Lyons, 2008). Of course, Lyons, 2008 should have 100% similarity for itself. The next highest green bar for this paper also happens to be Lyons, 2008. The similarity score is a bit lower (but sitll very high) because stopwords (e.g. the, a, an, them) have been removed.  
