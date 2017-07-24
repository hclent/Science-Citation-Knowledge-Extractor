from flask import Flask
import logging
from sqlalchemy import create_engine, MetaData, Table

app = Flask(__name__, static_url_path='/hclent/Webdev-for-bioNLP-lit-tool/flask/static')
app.config.from_pyfile('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/configscke.cfg', silent=False) #pass abs path

logging.basicConfig(filename='.app.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

def connect_db():
	engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'])
	return engine


def connection():
	conn = engine.connect()
	logging.info("NEW DB CONNECTION :')")
	return conn


def load_tables():
	metadata = MetaData(bind=engine) #init metadata. will be empty
	metadata.reflect(engine) #retrieve db info for metadata (tables, columns, types)
	inputPapers = Table('inputPapers', metadata)
	citations = Table('citations', metadata)
	queries = Table('queries', metadata)
	annotations = Table('annotations', metadata)
	return inputPapers, citations, queries, annotations


logging.info("DB INITIALIZING ... ")
engine = connect_db()
inputPapers, citations, queries, annotations = load_tables()
logging.info("DB INITIALIZEED !!! ")



