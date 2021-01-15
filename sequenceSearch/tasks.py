from __future__ import absolute_import, unicode_literals
from . import sequencetools
import os
from django.conf import settings
import sys

# Set clusterCAD settings for celery app
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'clusterCAD.settings')
sys.path.append('../')

# Import celery app 
from clusterCAD.celery import app

# This is a function wrapper for blast that is needed for Celery
# The actual function called in located in sequencetools.py
@app.task()
def blast(**kwargs):
	return sequencetools.blast(**kwargs)
   
