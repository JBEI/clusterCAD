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

#for corrfile in ['./data/corrections/modified/BGC0000055.1.json',
#                 './data/corrections/modified/BGC0000054.1.json',
#                 './data/corrections/modified/BGC0000416.1.json',
#                 './data/corrections/modified/BGC0000031.1.json',
#                 './data/corrections/modified/BGC0000430.1.json',
#                 ]:

#for corrfile in ['./data/corrections/modified/BGC0001394.1.json',
#                 './data/corrections/modified/BGC0001095.1.json',
#                 './data/corrections/modified/BGC0001024.1.json',]: # test with select clusters
for corrfile in correctionlist:
    
    acc = os.path.basename(corrfile).strip('.json')
    print('Correcting cluster %s.'  % acc)
    #try:
    cluster = pks.models.Cluster.objects.get(mibigAccession=acc)
    cluster.correctCluster(corrfile)
    #except MultipleObjectsReturned:
    #    # use some other method to identify object, see if database fixes all clusters, or skip (continue)
    #    continue

    # if this cluster has a corrections file we will consider it
    # 'reviewed' and set the boolean value accordingly
    cluster.reviewed = True
    cluster.save()
