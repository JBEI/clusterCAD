from __future__ import absolute_import, unicode_literals
from . import sequencetools
import os
from django.conf import settings
import sys

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'clusterCAD.settings')

sys.path.append('../')


from clusterCAD.celery import app


@app.task()
def blast(query, 
        evalue=10.0, 
        max_target_seqs=10, 
        sortOutput=True,
        database=os.path.join(
                settings.BASE_DIR, 
                'pipeline', 'data', 'blast', 'clustercad_subunits_reviewed',
            ),
    ):
    return sequencetools.blast(query, evalue, max_target_seqs, sortOutput, database)
   
