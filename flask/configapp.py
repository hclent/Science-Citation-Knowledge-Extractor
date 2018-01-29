from flask import Flask
import logging, os
from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy import exc, event, select

############## !!! IMPORTANT  !!!  IMPORTANT !!! ##################################################
'''
In order to run this code you **MUST EDIT THIS FILE WITH ABSOLUTE PATHS** for both static_url_path and 
the path to your config file. 
'''
app = Flask(__name__, static_url_path='/usr/src/Science-Citation-Knowledge-Extractor/flask/static') #pass abs path
app.config.from_pyfile('/usr/src/Science-Citation-Knowledge-Extractor/configscke.cfg', silent=False) #pass abs path
###################################################################################################


logging.basicConfig(filename='.app.log',level=logging.DEBUG)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

def connect_db():
	engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'], pool_size=20, pool_recycle=3600)
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

#pool_pre_ping does not work
@event.listens_for(engine, "engine_connect")
def ping_connection(connection, branch):
	if branch:
		logging.info(str(branch) + "is a branch. sub-connection of a connection. don't pre-ping")
		# "branch" refers to a sub-connection of a connection,
		# we don't want to bother pinging on these.
		return

	# turn off "close with result".  This flag is only used with
	# "connectionless" execution, otherwise will be False in any case
	save_should_close_with_result = connection.should_close_with_result
	connection.should_close_with_result = False

	try:
		# run a SELECT 1.   use a core select() so that
		# the SELECT of a scalar value without a table is
		# appropriately formatted for the backend
		connection.scalar(select([1]))
		logging.info("connection scalar: " + str(connection.scalar(select([1]))))
	except exc.DBAPIError as err:
		logging.info(err)
		# catch SQLAlchemy's DBAPIError, which is a wrapper
		# for the DBAPI's exception.  It includes a .connection_invalidated
		# attribute which specifies if this connection is a "disconnect"
		# condition, which is based on inspection of the original exception
		# by the dialect in use.
		if err.connection_invalidated:
			logging.info("error" + str(err.connection_invalidated))
			# run the same SELECT again - the connection will re-validate
			# itself and establish a new connection.  The disconnect detection
			# here also causes the whole connection pool to be invalidated
			# so that all stale connections are discarded.
			connection.scalar(select([1]))
		else:
			raise
	finally:
		# restore "close with result"
		connection.should_close_with_result = save_should_close_with_result
