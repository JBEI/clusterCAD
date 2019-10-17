from __future__ import absolute_import, unicode_literals
import os
from celery import Celery
from django.conf import settings

# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'clusterCAD.settings')

# Create celery app called ClusterCAD
app = Celery('ClusterCAD')

# Using a string here means the worker doesn't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')

# Load task modules from all registered Django app configs.
app.autodiscover_tasks()


# Debugging tasks
# @app.task
# def add11(x, y):
# 	return x+y
# @app.task()
# def debug_task(self):
#    print('Request: {0!r}'.format(self.request))
