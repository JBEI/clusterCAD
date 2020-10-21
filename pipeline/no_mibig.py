#!/usr/bin/python3

import glob
import os, sys
import json
import pickle

from Bio import SeqIO

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")


mibigpath = './data/mibig/raw' 
antismashpath = './data/antismash/split'


#for accession in ['BGC0001348']: # Debug with Borreledin, 31
#    clusterfile = os.path.join(antismashpath, accession + '.embl')
for file in glob.glob(os.path.join(antismashpath, '*.embl')):
    record = SeqIO.read(file, "embl")
    #print('\n'.join(dir(record)))
    #print(record.annotations)
    #print(record.features[-1])
    #TODO look only at cds annotations. look at sec_met and look for the domains, eg PKS_AT, PKS_KS, etc.
    #first identify all modules in pks and then pick the earliest
    #extension modules always have ks, at, acp, loading modules might not... We can look for at least an ACP, and then one search for something with all 3(ks, at,acp). The ones with 3 are probably extension, but if 
    #idea: if it has at least 1, then it's a pks module, if it has 3 it's a extension, so if it has <3 it's a loading...?


    print('processing', file, "...\n\n\n")
    for feature in record.features:
        #print(feature.qualifiers.values())
        if (feature.type == 'CDS' and 'sec_met' in feature.qualifiers):
            print('hello')
            #print((feature.qualifiers['sec_met']))
        #print(list(filter(lambda value: print(value), feature.qualifiers.values())))
        #print((feature.location))
        pass

