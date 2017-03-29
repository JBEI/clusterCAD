#!/usr/bin/python

import os, sys
import glob

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

correctionpath = './data/corrections'
correctiontargets = ['BGC0000029.1', 'BGC0000031.1', 'BGC0000042.1', 
                     'BGC0000093.1', 'BGC0000097.1', 'BGC0000115.1', 
                     'BGC0000165.1', 'BGC0001072.1', 'BGC0001381.1']

for acc in correctiontargets:
    cluster = pks.models.Cluster.objects.get(mibigAccession=acc)
    cluster.clusterJSON(correctionpath)
