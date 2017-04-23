#!/usr/bin/python3

import os
import sys
import pickle

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models
import compounddb.models

# Function to write pickled dictionaries containing sequence analysis information

def file_to_dict(dirname, filename, val_type='str'):
    assert val_type in ['str', 'int']
    # Initialize sequence dict
    sequence_dict = {}
    with open(os.path.join(dirname, filename), 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    key = lines[0][1:]
    sequence = []
    for line in lines[1:]:
        if line[0] == '>':
            sequence_dict[key] = sequence
            key = line[1:]
            sequence = []
        else:
            if val_type == 'str':
                sequence.extend(list(line))
            else:
                sequence.extend([int(v) for v in line.split()])
    return sequence_dict

# Directory containing sequence data
dirname = './data/aa_sequence_analysis'

acc_dict = file_to_dict(dirname, 'dataset.acc', 'str')
acc20_dict = file_to_dict(dirname, 'dataset.acc20', 'int')
ss_dict = file_to_dict(dirname, 'dataset.ss', 'str')
ss8_dict = file_to_dict(dirname, 'dataset.ss8', 'str')

pickle.dump(acc_dict, open(os.path.join(dirname, 'acc_dict.p'), 'wb'))
pickle.dump(acc20_dict, open(os.path.join(dirname, 'acc20_dict.p'), 'wb'))
pickle.dump(ss_dict, open(os.path.join(dirname, 'ss_dict.p'), 'wb'))
pickle.dump(ss8_dict, open(os.path.join(dirname, 'ss8_dict.p'), 'wb'))

# Insert sequence information into database
for cluster in pks.models.Cluster.objects.all():
    for subunit in pks.cluster.subunits():
        # Reconstruct reference string
        reference = subunit.genbankAccession
        reference += '_'
        name = subunit.name.split('_')
        reference += name[0]
        if name[-1] != name[0]:
            reference += name[-1]
            order = subunit.order
            if order != 0:
                reference += '_%s' %(order)
        # Add data to database
        try:
            subunit.acc = acc_dict[reference]
            subunit.acc20 = acc20_dict[reference]
            subunit.ss = ss_dict[reference]
            subunit.ss8 = ss8_dict[reference]
        except Exception:
            raise Exception('Missing data!')
