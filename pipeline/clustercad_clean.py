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
    except Exception as e:
        print(e)
        print('FAILED COMPUTE PRODUCT...', cluster)#, cluster.subunits()[0].modules()[0].domains())
        for subunit in cluster.subunits():
            print(subunit.modules())
            for m in subunit.modules():
                print(m.domains())


        cluster.delete()
        print('Deleted %s: %s because no computed product' %(cluster.mibigAccession, cluster.description))
        continue
    if nmodules < 3:
        print('deleted %s: %s becuase less than 3 modules' %(cluster.mibigAccession, cluster.description))
        cluster.delete()

# Delete clusters with no computable product
for cluster in pks.models.Cluster.objects.all():
    try:
        cluster.computeProduct()
    except:
        print('Deleted %s: %s because no computed product... again.' %(cluster.mibigAccession, cluster.description))
        cluster.delete()

