from celery import task

@task
def greeting():
    print("hello i have no idea what i'm doing")