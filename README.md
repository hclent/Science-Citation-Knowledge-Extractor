# Science Citation Knowledge Extractor (SCKE)


## About SCKE
* SCKE is an open source tool that helps biomedical researchers understand how their work is being used by others, by analyzing the content in papers that cite them. This tool uses natural language processing and machine learning to extract the prominent themes and concepts discussed in the citing documents. By seeing what kinds of topics are discussed by citing articles, researchers can better understand how their work influences their peers and various disciplines in science. Additionally, SCKE allows biomedical researchers to explore other statistics about the publications that site them, such as where citations are published (Journals), the distribution of keywords (Keywords), the similarity of papers to each other (Clustering), the similarity of papers to other famous works (TextCompare), and general statistics about the citations (Statistics).
* Built with `Python3` and `Flask` using `Biopython`, `py-processors`, `NLTK`, `Scikit-Learn`, `Numpy`, `Gensim`, `Fasttext`, `Plotly`, `D3`, and more.
* SCKE proudly leverages the BioNLP Processor created by the Computational Language Understanding Lab. See [here](https://github.com/clulab/processors) for original code. See [here](https://github.com/myedibleenso/py-processors) for python wrappper library, `py-processors`.

## Instructions for using the site
* Video to come :)

## Instructions for setting up the site

**Requirements:**

* Python3 ([Anaconda](https://www.continuum.io/))
* [MySQL server](https://www.mysql.com/)

**Step 1: Retrieve the source code**

1. `git clone` to run this repository

**Step 2: Setting up Python**

In order to set up the site, you will need the same Python environment that the site uses to run. Here are instructions for how to re-create this environment.

1. Create a Python environment with the repository's `requirements.txt` libraries (e.g `conda create -n scke python=3.4 anaconda --file requirements.txt`)
2. Activate the environment (`source activate scke`)

**Step 3: Creating a config file**

Next you will need to create a `configscke.cfg` file, that tells the files in the app where to look for the database and cached data. For complete instructions, please see the example `configExample.cfg`.

**Step 4: Setting up MySQL database**

1. Create a new database for scke data
2. See table configuration notes in [database_management.py](https://github.com/hclent/Webdev-for-bioNLP-lit-tool/blob/master/flask/database_management.py#L583)

**Step 5: Add your credentials for Entrez**

Entrez is the API for accessing PubMed(Central). We use Entrez as part of the Biopython package to scrape the citations. You should always tell Entrez who you are.

1. Simply edit [this line of Entrez_IR.py](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/Entrez_IR.py#L17) with your email.

**Step 6: Other files you will need**

These files live on CyVerse's Discovery environment (link to come)

1. The py-processors NLP Server jar
    * This jar is used to annotate the texts. It is initialized to use 3G of memory but you can increase this in the `configscke.cfg` file.
    * If you you are running SCKE with either large documents or on a large amount of documents 3G will not be enough memory for the jar file! We give SCKE 100G of memory to run it in production. If you experience problems, try bumpting the memory up to 5G, 10G, or 25G.
2. Our trained `FastText` word vectors `.vec` file

**Step 7: Run SCKE**

Once you have cloned the repo, made your new Python environment, created a config file, and set up the MySQL database, you are almost ready to run SCKE!

1. Activate your conda environment (`source activate scke`)
2. Modify [this line of app.py](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/app.py#L1110).
    * Comment out the line with `run_simple`
        * `run_simple` is for running the app with Apache and uwsgi
    * Uncomment the line below it, with `app.run()`
3. Navigate to the `flask` directory of the repository (`cd flask`)
4. Run `python app.py`

**Debugging?**

If you have any problems, feel free to contact me! But first, take a look instead of your `.app.log` file to see where things might have gone amiss.