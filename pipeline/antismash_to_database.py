#!/usr/bin/python3

import os, sys
import json
import pickle

from Bio import SeqIO
from antismash_to_database_functions import mibigSubtypes, filterModularTypeI, getNRPS, enterCluster

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

##########################################################
# Identifying valid type I modular PKSs and NRPS from MIBiG files #
##########################################################

print('Analyzing contents of MIBiG database.')
mibigpath = './data/mibig/raw' 
antismashpath = './data/antismash/split'
print('Set of PKS subtypes found in MIBiG: %s' %(mibigSubtypes(mibigpath)))
mibigsubtypes = set(['Modular type I', 'Modular Type I', 'Type I'])
print('Set of PKS subtypes recognized for inclusion in ClusterCAD: %s' %(mibigsubtypes))
t1pksmibigaccessions = filterModularTypeI(mibigpath, mibigsubtypes)
mibigaccessions = getNRPS(mibigpath, t1pksmibigaccessions)
print('\n')

##################################
# ClusterCAD database generation #
##################################

print('Resetting ClusterCAD database.')
[cluster.delete() for cluster in pks.models.Cluster.objects.all()]
print('ClusterCAD database reset.')

# Assumes that chemical structures have already been aggregated
allknowncompounds = pickle.load(open('./data/compounds/all_known_products.p', 'rb'))

#for accession in ['BGC0000031']: # Debug with Borreledin
#for accession in ['BGC0001394','BGC0001024','BGC0000416', 'BGC0001095', 'BGC0000430', 'BGC0000367']: # Debug NRPS pipeline with Phenalamide and friends
for accession in mibigaccessions:
    # Use accession number to get paths to MIBiG and antiSMASH files
    mibigfile = os.path.join(mibigpath, accession + '.json')
    clusterfile = os.path.join(antismashpath, accession + '.embl')

    # Read antiSMASH annotations for cluster
    try:
        record = SeqIO.read(clusterfile, "embl")
        description = record.description.replace(' biosynthetic gene cluster', '')
    # If file is missing, we skip the cluster
    except FileNotFoundError:
        print('Missing file: %s' %clusterfile)
        continue

    # Get compound information
    try:
        compound = allknowncompounds[accession]
    # If compound is missing, we skip the cluster
    except KeyError:
        print('Missing compound %s: %s.' %(accession, description))
        continue
    knownproductsmiles = compound[0][0]
    knownproductsource = compound[1]

    # Enter information in ClusterCAD database
    try:
        cluster = pks.models.Cluster(
            genbankAccession=record.annotations['comment'].split()[-1].strip().strip('.'), \
            mibigAccession=record.id, \
            description=description, \
            sequence=record.seq,
            knownProductSmiles=knownproductsmiles,
            knownProductSource=knownproductsource
            )
        print(record)
        cluster.save()
        # Processes subunits and modules belonging to cluster
        enterCluster(cluster, record, mibigfile)
        print('Processed cluster %s: %s.' %(record.id, record.description))
        cluster.computeProduct()
        print('Pregenerated cluster products.')
    except Exception as e:
        print(e)
        pass
