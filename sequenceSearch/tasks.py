from __future__ import absolute_import, unicode_literals
from . import sequencetools
import os
from django.conf import settings
import sys
from django.core.cache import cache

# Set clusterCAD settings for celery app
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'clusterCAD.settings')
sys.path.append('../')

# Import celery app 
from clusterCAD.celery import app

# This is a function wrapper for blast that is needed for Celery
# The actual function called in located in sequencetools.py
@app.task()
def blast(**kwargs):
    cachekey = tuple(sorted(kwargs.items()))
    alignments = cache.get(cachekey)
    if alignments is None:
        alignments = sequencetools.blast(**kwargs)
        cache.set(cachekey, alignments, 60 * 60 * 24 * 31) # cache for 31 days

    return alignments

   
