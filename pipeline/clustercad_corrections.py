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

for corrfile in correctionlist:
    
    acc = os.path.basename(corrfile).strip('.json')
    print('Correcting cluster %s.'  % acc)
    cluster = pks.models.Cluster.objects.get(mibigAccession=acc)
    cluster.correctCluster(corrfile)

    # if this cluster has a corrections file we will consider it
    # 'reviewed' and set the boolean value accordingly
    cluster.reviewed = True
    cluster.save()

print("\nAll cluster corrections applied. ")
