# Science Citation Knowledge Extractor (SCKE)


## About SCKE
* SCKE is an open source tool that helps biomedical researchers understand how their work is being used by others, by analyzing the content in papers that cite them. This tool uses natural language processing and machine learning to extract the prominent themes and concepts discussed in the citing documents. By seeing what kinds of topics are discussed by citing articles, researchers can better understand how their work influences their peers and various disciplines in science. Additionally, SCKE allows biomedical researchers to explore other statistics about the publications that site them, such as where citations are published (Journals), the distribution of keywords (Keywords), the similarity of papers to each other (Clustering), the similarity of papers to other famous works (TextCompare), and general statistics about the citations (Statistics).
* Built with `Python3` and `Flask` using `Biopython`, `py-processors`, `NLTK`, `Scikit-Learn`, `Numpy`, `Gensim`, `Fasttext`, `Plotly`, `D3`, and more.
* SCKE proudly leverages the BioNLP Processor created by the Computational Language Understanding Lab. See [here](https://github.com/clulab/processors) for original code. See [here](https://github.com/myedibleenso/py-processors) for python wrappper library, `py-processors`.

## Instructions for using the site
* Video to come :)

## Docker container for running the site
* Docker container(s) to come!


## Instructions for setting up the site without a Docker container

**Requirements:**

* Python3 ([Anaconda](https://www.continuum.io/))
* [MySQL](https://www.mysql.com/)
* Java
* [Bower](https://bower.io/)

**Step 1: Install MySql and configure database**

1. Make sure MySql is installed. (We use version `5.6.33`)
2. Create a new database for scke data with `create database scke`
2. Create the necessary tables. See table configuration notes in [database_management.py](https://github.com/hclent/Webdev-for-bioNLP-lit-tool/blob/master/flask/database_management.py#L583) for instructions.

**Step 2: Install Java**

Make sure you have Java installed!

**Step 3: Retrieve the source code**

Clone (or fork!) the project and move into the main directory

1. `git clone
https://github.com/hclent/Science-Citation-Knowledge-Extractor.git`
2. `cd Science-Citation-Knowledge-Extractor`.

**Step 4: Get the javascript with Bower**

A lot of the Javascript is already provided in the repo, but there are still a few packages you will need. `bower install bower.json` to get them.

**Step 5: Download the supplementary files (You need these!)**

These supplementary files are stored on CyVerse's [Discovery Environment](http://www.cyverse.org/discovery-environment)

1. Install [iCommands](https://pods.iplantcollaborative.org/wiki/display/DS/Using+iCommands)
2. (more instructions to come  ....)

The files in the Discovery Environment that you will need are:

* The py-processors  BioNLP Server jar, `processors-server.jar`
    * This jar is used to annotate the texts. It is initialized to use 10G of memory but you can increase by editing [this line of multi_preprocess.py](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/multi_preprocess.py#L30).
    * If you you are running SCKE with either large documents or on a large amount of documents 10G will not be enough memory for the jar file! We give SCKE 100G of memory to run it in production. If you experience problems, try bumpting the memory up to 25G, 50G, or more.
* Our trained `FastText` word vectors file, `scke_lower.vec`

Appart from the Discovery Environment, we also make use of `nltk`'s stopwords. You can download these easily by entering the python shell and

```
import nltk
nltk.download()
```
A GUI will come up. Type `d stopwords`, download the stopwords, then `q` to quit and `exit()`.

**Step 6: Setting up Python**

In order to set up the site, you will need the same Python environment that the site uses to run. Here are instructions for how to re-create our environment.

1. Install ([Anaconda](https://www.continuum.io/)) for Python3
2. Create a Python environment with the repository's `requirements.txt` libraries (e.g `conda create -n scke python=3.4 anaconda --file conda_requirements.txt`)
3. Activate the environment (`source activate scke`)
4. Install the rest of the requirements in the environment with `pip install -r pip_requirements.txt`

**Step 7: Creating and loading the config file**

1. Next you will need to create a config file, that tells the files in the app where to look for the database and cached data, etc. For complete instructions, please see the example `configExample.cfg`.
2. Make sure that [this line of](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/configapp.py#L7)`configapp.py` is pointing to the correct location of your config file!
    * NB: Also while you're updating this, [the line above it](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/configapp.py#L6) should also be pointing to the `static` directory.

**Step 8: Add your credentials for Entrez**

Entrez is the API for accessing PubMed(Central). We use Entrez as part of the Biopython package to scrape the citations. You should always tell Entrez who you are.

1. Simply edit [this line of Entrez_IR.py](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/Entrez_IR.py#L17) with your email.


**Step 9: Run SCKE**

Once you have accomplished steps 1-8, you are almost ready to run SCKE!

1. Make sure your conda environment is active (`source activate scke`)
2. Modify [this line of app.py](https://github.com/hclent/Science-Citation-Knowledge-Extractor/blob/master/flask/app.py#L1110).
    * Comment out the line with `run_simple`
        * `run_simple` is for running the app with Apache and uwsgi
    * Uncomment the line below it, with `app.run()`
3. Navigate to the `flask` directory of the repository (`cd flask`)
4. Run `python app.py`

It takes a few minutes to start up since it must establish a connection with our `processors-server.jar` and load the `scke_lower.vec` file. You can `tail -F .app.py` to watch the log file go.

To view the website go to either `0.0.0.0:5000\home` if you are running it locally or if you are running this on a vm go to `<THE_IP_ADDRESS>:5000/home`

**Debugging?**

If you have any problems, feel free to contact me! But first, take a look instead of your `.app.log` file to see where things might have gone amiss.