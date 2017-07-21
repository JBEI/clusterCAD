#!/usr/bin/python3

import os, sys
import glob

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

correctionpath = './data/corrections/modified'
correctionlist = glob.glob(os.path.join(correctionpath, '*.json'))

for corrfile in correctionlist:
    acc = os.path.basename(corrfile).strip('.json')
    print('Correcting cluster %s.'  % acc)
    cluster = pks.models.Cluster.objects.get(mibigAccession=acc)
    cluster.correctCluster(corrfile)
