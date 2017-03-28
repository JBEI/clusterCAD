#!/usr/bin/python

import os, sys

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

# Delete clusters with less than three modules
for cluster in pks.models.Cluster.objects.all():
    nmodules = 0
    for subunit in cluster.subunits():
        nmodules += len(subunit.modules())
    if nmodules < 3:
        print('%s: %s' %(cluster.mibigAccession, cluster.description))
        cluster.delete()

# Delete clusters with no computable product
for cluster in pks.models.Cluster.objects.all():
    try:
        cluster.computeProduct()
    except:
        print('%s: %s' %(cluster.mibigAccession, cluster.description))
        cluster.delete()
