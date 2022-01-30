#!/usr/bin/python3

import os, sys
import glob

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

correctionpath = './data/corrections'
# correctiontargets = ['BGC0000029.1', 'BGC0000031.1', 'BGC0000042.1', 
#                     'BGC0000093.1', 'BGC0000097.1', 'BGC0000165.1', 
#                     'BGC0001381.1']

correctiontargets = ['BGC0001296.1']
#correctiontargets = [x.mibigAccession for x in pks.models.Cluster.objects.all()]

for acc in correctiontargets:
    print('Getting template for cluster %s' %(acc))
    cluster = pks.models.Cluster.objects.get(mibigAccession=acc)
    cluster.clusterJSON(correctionpath)
