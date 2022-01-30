#!/usr/bin/python3

import os, sys
import pandas as pd

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models


### new pipeline already deletes subunits without modules, clusters with less than 3 modules OK as unreviewed ###
'''
# Delete clusters with less than three modules
for cluster in pks.models.Cluster.objects.all():
    nmodules = 0
    for subunit in cluster.subunits():
        nsubmodules = len(subunit.modules())
        if nsubmodules == 0:
            subunit.delete()
        nmodules += nsubmodules
    # Recompute product once invalid subunits have been deleted
    try:
        cluster.computeProduct(recompute=True)
    except:
        print(cluster.architecture())
        cluster.delete()
        print('could not recompute product, deleted %s: %s' %(cluster.mibigAccession, cluster.description))
        continue
    #if nmodules < 3:
    #    print('less than 3 modules, deleted %s: %s' %(cluster.mibigAccession, cluster.description))
    #    cluster.delete()
'''
# Delete clusters with no computable product
for cluster in pks.models.Cluster.objects.all():
    try:
        cluster.computeProduct()
    except:
        print('could not compute product, deleted %s: %s' %(cluster.mibigAccession, cluster.description))
        cluster.delete()

