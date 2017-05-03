from datetime import timedelta

CELERY_IMPORTS = ('task')
CELERY_IGNORE_RESULT = False
BROKER_HOST = '127.0.0.1'
BROKER_PORT = 5672
BROKER_URL = 'amqp://'

CELERYBEAT_SCHEDULE = {
    'sayhello': {
        'task': 'task.greeting',
        'schedule': timedelta(seconds=10),
    },
}