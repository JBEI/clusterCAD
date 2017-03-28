#!/usr/bin/python

import os, sys
import json
import pickle

from Bio import SeqIO

from antismash_to_database_functions import mibigSubtypes, filterModularTypeI, enterCluster
sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models
import compounddb.models

##########################################################
# Identifying valid type I modular PKSs from MIBiG files #
##########################################################

print('Analyzing contents of MIBiG database.')
mibigpath = './data/mibig/raw' 

antismashpath = '../notebooks/mibig/antismash/split_files'





#antismashpath = './data/antismash/split'
print('Set of PKS subtypes found in MIBiG: %s' %(mibigSubtypes(mibigpath)))
mibigsubtypes = set(['Modular type I', 'Modular Type I', 'Type I'])
print('Set of PKS subtypes recognized for inclusion in ClusterCAD: %s' %(mibigsubtypes))
mibigaccessions = filterModularTypeI(mibigpath, mibigsubtypes)
print('\n')

##################################
# ClusterCAD database generation #
##################################

print('Resetting ClusterCAD database.')
[cluster.delete() for cluster in pks.models.Cluster.objects.all()]
print('ClusterCAD database reset.')

# Assumes that chemical structures have already been aggregated
allknowncompounds = pickle.load(open('./data/compounds/all_known_products.p', 'rb'))

for accession in mibigaccessions:
    # Use accession number to get paths to MIBiG and antiSMASH files
    mibigfile = os.path.join(mibigpath, accession + '.json')
    clusterfile = os.path.join(antismashpath, accession + '.embl')

    # Read antiSMASH annotations for cluster
    record = SeqIO.read(clusterfile, "embl")

    # Get compound information
    try:
        compound = allknowncompounds[record.id.split('.')[0]]
    # If compound is missing, we skip the cluster
    except KeyError:
        continue
    knownproductsmiles = compound[0][0]
    knownproductsource = compound[1]

    # Enter information in ClusterCAD database
    try:
        cluster = pks.models.Cluster(
            genbankAccession=record.annotations['comment'].split()[-1].strip().strip('.'), \
            mibigAccession=record.id, \
            description=record.description.replace(' biosynthetic gene cluster', ''), \
            sequence=record.seq,
            knownProductSmiles=knownproductsmiles,
            knownProductSource=knownproductsource
            )
        cluster.save()
        # Processes subunits and modules belonging to cluster
        enterCluster(clusterfile, clusterfile, mibigfile)
    except Exception as e:
        pass
