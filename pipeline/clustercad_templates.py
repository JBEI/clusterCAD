#!/usr/bin/python3

import os, sys
import glob

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

correctionpath = './data/corrections'

#correctiontargets = ['BGC0001792.1']
correctiontargets = [x.mibigAccession for x in pks.models.Cluster.objects.all()]

for acc in correctiontargets:
    print('Getting template for cluster %s' %(acc))
    cluster = pks.models.Cluster.objects.get(mibigAccession=acc)
    cluster.clusterJSON(correctionpath)
