#!/usr/bin/python3

import glob
import re
import os, sys
import json
import pickle

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

import django
django.setup()
import pks.models
from Bio import SeqIO

mibigpath = './data/mibig/raw' 
antismashpath = './data/antismash/split'
newantismashpath = './pipeline/antismash_db_sample'

for file in glob.glob(os.path.join(newantismashpath, '*.gbk')):
    record = SeqIO.read(file, "genbank")


    print('processing', file, "...\n")
    for feature in record.features:
        """
        feature.qualifiers is a list with format as follows:
        [
            'Type: xxxxx',
            'Domains detected: domainA (E-value:123, bitscore: 123, seeds: 123); domainB (E-value...'
        ]
        """
        if (feature.type != 'CDS' or 'sec_met' not in feature.qualifiers):
            continue
        if ('gene' not in feature.qualifiers.keys() and 'product' not in feature.qualifiers.keys()):
            continue
        print(feature.qualifiers)
        continue
        if (not sec_met_type): 
            # if no regex match, it probably looks something like this:
            # NRPS/PKS subtype: PKS/NRPS-like protein 
            print('skipped', feature.qualifiers['protein_id'])
            continue
        if sec_met_type[1] != 't1pks':
            #check if the xxx in Type: xxx is 't1pks'.
            #continue
            pass
        #at this stage, we now have a record feature that is labelled t1pks by antismash.
        """sec_met_domains = ' '.join(feature.qualifiers['sec_met'][1:])
        contains_ks = 'PKS_KS' in sec_met_domains # there is also a mod_KS but does that count?
        contains_at = 'PKS_AT' in sec_met_domains
        contains_acp = 'ACP' in sec_met_domains
        #print(("protein_id" in feature.qualifiers) and feature.qualifiers['protein_id'], contains_ks, contains_at, contains_acp)"""

        if 'gene' in feature.qualifiers.keys(): 
            genename = feature.qualifiers['gene'][0]
        elif 'product' in feature.qualifiers.keys():
            genename = feature.qualifiers['product'][0]
        subunit = pks.models.Subunit(cluster=TODO,
                                     genbankAccession=feature.qualifiers['protein_id'][0],
                                     name=genename,
                                     start=feature.location.start.position + 1,
                                     stop=feature.location.end.position
                                     sequence=feature.qualifiers['translation'][0])
        subunit.save()

