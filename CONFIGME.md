##Setting up `configscke.cfg`

Here we will add the full paths for important files and directores that SCKE uses.

First create the config file `touch configscke.cfg` file.
Then `vi configscke.cfg` to add your configuration settings

* `SECRET_KEY = 'your secret key'`
* `SQLALCHEMY_DATABASE_URL = 'mysql:://username@password@localhost:port/db_name?charset=utf8mb4'`
    * the `?charset=utf8mb4` ensures that there are no issues with unicode in the database
* `SCKE_URL_RESULTS = 'https://www.yourwebsite.com/results'`
* `PATH_TO_LOG = '/the/path/to/your/repo/flask/.app.log'`
    * The log is used for debugging and is located in the project's `flask` directory
* `PATH_TO_JAR = 'your/home/anacond3/envs/scke/py34/lib/python3.4/site-packages/processors/processor-server.jar'`
     * This points to a file living inside your Anaconda distribution of Python from the `conda env`. This `processor-server.jar` is necessary for annotating our `.txt` documents into "BioDocs", that contain information about dependency parsing, biomedical named entity recognition, and more.
* `PATH_TO_CACHE` is where the scraped `.txt` documents, the annotated biodocs (`.json`), and the special query-related `.pickle` files are located
* `PATH_TO_STATIC` is the path to `static` inside the `flask` dir. This is a typical file to have for any `flask` project.
* `PATH_TO_FRAPHS` and `PATH_TO_CLUSTERM` should point to folders in `static` where the cached files for the force directed graphs and clustermaps are located, respectively.
* `PATH_TO_FASTTEXT_MODEL` should point to a trained `.vec` file
* `PATH_TO_JOURNALS`, `PATH_TO_LDA`, `PATH_TO_LSA` should point to where you intend to store these data visualizations


When are you are done, hit `esc` and then `:wq` to save and quit