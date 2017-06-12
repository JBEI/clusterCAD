#!/usr/bin/python3

import os, sys
import pandas as pd

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

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
        cluster.delete()
        print('%s: %s' %(cluster.mibigAccession, cluster.description))
        continue
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

# read in clusters we decided not to fix
rejectclusters = pd.read_csv("data/corrections/modified/drop_clusters.txt", comment="#", header=None)[0].tolist()
rejectclusters = [x.strip() for x in rejectclusters]

# Delete clusters we've decided not to fix
for acc in rejectclusters:
    cluster = pks.models.Cluster.objects.get(mibigAccession=acc)
    print('%s: %s' %(cluster.mibigAccession, cluster.description))
    cluster.delete()
