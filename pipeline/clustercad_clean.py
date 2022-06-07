#!/usr/bin/python3

import os, sys
import pandas as pd

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models


# Empty subunits are deleted by the pipeline as they are generated but clean checks again
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
        print('Could not recompute product, deleted %s: %s' %(cluster.mibigAccession, cluster.description))
        continue
    # Delete clusters with no modules
    if nmodules == 0:
        print('No modules detected, deleted %s: %s' %(cluster.mibigAccession, cluster.description))
        cluster.delete()

# Delete clusters with no computable product
for cluster in pks.models.Cluster.objects.all():
    try:
        cluster.computeProduct()
    except:
        print('Could not compute product, deleted %s: %s' %(cluster.mibigAccession, cluster.description))
        cluster.delete()

print("\nDatabase cleanup completed. ")