# Science Citation Knowledge Extractor (SCKE)
* UNDER ACTIVE CONSTRUCTION
* COMING SOON!

## BASICS
* This will be the web interface for the bioNLP lit tool
* Built with Python3 and Flask

## Instructions for using the site
* Video to come :)

## Instructions for setting up the site

**Requirements:**

* Python3 ([Anaconda](https://www.continuum.io/))
* [MySQL server](https://www.mysql.com/)
* [Apache Web Sever](https://en.wikipedia.org/wiki/Apache_HTTP_Server)
    * Only necessary if you want to run it on a specific website domain.
    * If you want to run the app on your `localhost`, then Apache is not neeeded.

**Step 1: Retrieve the source code**

1. `git clone` or `fork` this repository

**Step 2: Setting up Python**

In order to set up the site, you will need the same Python environment that the site uses to run. Here are instructions for how to re-create this environment.

1. Create a Python environment with the repository's `requirements.txt` libraries (e.g `conda create -n scke python=3.4 anaconda --file requirements.txt`)
2. Activate the environment (`source activate scke`)

**Step 3: Creating a config file**

Next you will need to create a `configscke.cfg` file, that tells the files in the app where to look for the database and cached data. For complete instructions, please see `CONFIGME.md`

**Step 4: Setting up MySQL database**

1. Create a new database for scke data
2. See table configuration notes in [database_management.py](https://github.com/hclent/Webdev-for-bioNLP-lit-tool/blob/master/flask/database_management.py#L578)

**Step 5: Run SCKE**
Once you have set up Apache, made your new Python environment, created a config file, and set up the MySQL database, you are ready to run SCKE!

1. Activate your conda environment (`source activate scke`)
2. Navigate to the `flask` directory of the repository (`cd flask`)

**If you want to run this on a specific website**

Run SCKE with the command: `uwsgi --socket YOUR_HOST:YOUR_PORT --protocol=https --manage-script-name --mount /test=/path/to/this/repo/flask/app.py --callable app -p 20`
    * the `-p 20` is there to give the app extra workers
    * an example of `YOUR_HOST:YOUR_PORT` could be `0.0.0.0:8888`


**If you want to run this on localhost:**
For running it on your localhost, you will need to:

1. Modify [this line of app.py](https://github.com/hclent/Webdev-for-bioNLP-lit-tool/blob/master/flask/app.py#L1052).
    * Comment the line with `run_simple`
    * Uncomment the line below it, with `app.run()`
2. Run `python app.py` 