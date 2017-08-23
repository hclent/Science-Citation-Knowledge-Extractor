from flask import Flask
from celery import Celery
import logging
from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy import exc, event, select

app = Flask(__name__, static_url_path='/hclent/Webdev-for-bioNLP-lit-tool/flask/static')
app.config.from_pyfile('/home/hclent/repos/Webdev-for-bioNLP-lit-tool/configscke.cfg', silent=False) #pass abs path

# TODO: make sure celery is configured correctly
# celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
# celery.conf.update(app.config)

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
	except exc.DBAPIError as err:
		# catch SQLAlchemy's DBAPIError, which is a wrapper
		# for the DBAPI's exception.  It includes a .connection_invalidated
		# attribute which specifies if this connection is a "disconnect"
		# condition, which is based on inspection of the original exception
		# by the dialect in use.
		if err.connection_invalidated:
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


#TODO: Celery tasks for asynch jobs, modal message updates, & sending e-mails
# @celery.task
# def greeting():
# 	print("hello i have no idea what i'm doing")
#
# task = greeting.apply_async(countdown=10)

# @celery.task(bind=True)
# def long_task(self):
# 	"""Background task that runs a long function with progress reports."""
# 	verb = ['Starting up', 'Booting', 'Repairing', 'Loading', 'Checking']
# 	adjective = ['master', 'radiant', 'silent', 'harmonic', 'fast']
# 	noun = ['solar array', 'particle reshaper', 'cosmic ray', 'orbiter', 'bit']
# 	message = ''
# 	total = random.randint(10, 50)
# 	for i in range(total):
# 		if not message or random.random() < 0.25:
# 			message = '{0} {1} {2}...'.format(random.choice(verb),
# 											  random.choice(adjective),
# 											  random.choice(noun))
# 		self.update_state(state='PROGRESS',
# 						  meta={'current': i, 'total': total,
# 								'status': message})
# 		time.sleep(1)
# 	return {'current': 100, 'total': 100, 'status': 'Task completed!',
# 			'result': 42}