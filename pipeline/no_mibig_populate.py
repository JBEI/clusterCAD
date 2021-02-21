#!/usr/bin/python3

"""
Populates from a directory specified in the first line of code.
The directory should contain gbk files from AntiSMASH 3.
This assumes that all files are already identified by AntiSMASH to be
Modular Type I PKS.
"""
ANTISMASH_FILES_DIR = '../antismash_db/'

import os, sys
import json
import pickle
import glob

from Bio import SeqIO
import antismash_to_database_functions as antismash_funcs
from antismash_to_database_functions import filterModularTypeI, enterCluster, processClusterSeqRecord, processSubunitModules, allowed_domains

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

filelist = glob.glob(os.path.join(ANTISMASH_FILES_DIR, '*.gbk'))

print(sys.argv)
#save_to_db if not a dry run
save_to_db = not ('-x' in sys.argv or '--no-save' in sys.argv or '--dry-run' in sys.argv)
clear_db = '-c' in sys.argv or '--clear' in sys.argv
debug = '-d' in sys.argv or '--debug' in sys.argv
accept_in_trans_modules = '--no-in-trans' not in sys.argv

def read_subunit(record):
    """Attempts to read the subunit information from
    a Cluster/gene record. Closely mimicks the processClusterSeqRecord
    function from the antismash db functions"""
    subunit_list = []
    for feature in record.features:
        if feature.type != 'CDS' or 'protein_id' not in feature.qualifiers.keys():
            continue
        subunit_info = {
            'protein_id': feature.qualifiers['protein_id'],
            'description': None,
            'genename': None,
            'translation': None,
            'location': None,
            'locus_tag': None,
        }
        if 'gene' not in feature.qualifiers.keys() and 'product' not in feature.qualifiers.keys():
            print('skipped subunit with no gene or product info.\n')
            continue
        # General information about the gene
        if 'product' in feature.qualifiers.keys():
            subunit_info['description'] = feature.qualifiers['product'][0]
        if 'gene' in feature.qualifiers.keys():
            subunit_info['genename'] = feature.qualifiers['gene'][0]
        elif 'product' in feature.qualifiers.keys():
            subunit_info['genename']= feature.qualifiers['product'][0]
        if 'sec_met' not in feature.qualifiers.keys():
            continue
        if 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag']:
            subunit_info['locus_tag'] = feature.qualifiers['locus_tag'][0]
        if len(feature.qualifiers['sec_met']) < 4:
            #print('\tsec_met field not big enough:\t', feature.qualifiers['sec_met'])
            continue
            """
            Information about domains is not completely and accurate captured in 
            the second entry, where it says "domains detected". We need to look 
            for "PKS/NRPS Domain:" entries which only happen at index greater than 4.
            """
        pks_subtype = feature.qualifiers['sec_met'][3]
        if pks_subtype in \
            ['NRPS/PKS subtype: Type I Modular PKS',
           'NRPS/PKS subtype: PKS-like protein',
           'NRPS/PKS subtype: PKS/NRPS-like protein',
           'NRPS/PKS subtype: Hybrid PKS-NRPS',
             ]:
            #print(feature.qualifiers['sec_met'])
            subunit_info['modules'] = antismash_funcs.processSubunitModules(feature.qualifiers['sec_met'])
            if debug:
                print('\t\t',feature.qualifiers['sec_met'])
                for order, module in subunit_info['modules'].items():
                    printable = ','.join(module.keys())
                    print(f'\t\t{order}: {printable}')
                print('-' * 200)
        elif pks_subtype == 'NRPS/PKS subtype: Type I Iterative PKS': 
            # Also process iterative PKS, with a note in the description
            subunit_info['description'] = 'Iterative PKS: ' + subunit_info['description']
            subunit_info['modules'] = antismash_funcs.processSubunitModules(feature.qualifiers['sec_met'], accept_in_trans_modules=True)
        else:
            print(f'annotation type "{pks_subtype}" not recognized and skipped.')
            continue
            #standalone information below not supported.
            #print('sec_met type not recognized, setting as satandalone')
            #print(feature.qualifiers['sec_met'])
            #subunit_info['modules'] = None
        subunit_info['location'] = feature.location.start.position + 1, feature.location.end.position #a tuple
        subunit_info['translation'] = feature.qualifiers['translation']
        subunit_list.append(subunit_info)
    return subunit_list

def process_subunits(subunits, cluster):
    """
    Given a list of subunits (Represented as dictionaries) returned from
    read_subunits, process them and save them to the actual database"""
    seen_genenames = set()
    counter = 1
    # first one uses this and then sets it false for later modules
    # TODO we have no idea if this is actually the loading module...
    is_loading_module = True
    for subunit in subunits:
        if not subunit['modules']: 
            # No modules, could be standalong, but not supported for now.
            continue
        # TAddress gene name duplicates by appendign a numebr at the end
        if subunit['genename'] in seen_genenames:
            subunit['genename'] = subunit['genename'] +  '_' + str(counter)
            counter += 1
        db_subunit_entry = pks.models.Subunit(cluster=cluster,
                                  genbankAccession=subunit['protein_id'][0], # TODO
                                  name = subunit['genename'] + ': ' + subunit['locus_tag'],
                                  start=subunit['location'][0],
                                  stop=subunit['location'][1],
                                  sequence=subunit['translation'],
                                  order=-1
                                  #TODO maybe needs ordering = false
                                  )
        if save_to_db:
            db_subunit_entry.save()
        print(f'\tsaved subunit entry {subunit["genename"]}')
        modules = subunit['modules']
        """

        module is the return value of processSubunitModules:
        keys: module number
        values: OrderedDict of domains in module. Within the OrderedDict, 
                 the key is the domain name, while the value is a length one or 
                 two list where the first element is a dictionary 
                 {start: value, stop: value} and the optional second element is 
                 a specificity dict if applicable
        example:
        {
            0: { #first module at location 0
                'KS': [ {'start': 0, 'stop': 1} ]
                'AT': [ {'start': 2, 'stop': 3} ]
            },
            1: { # a second module
                ...
            }
        }
        """
        for module_index, module_domains in modules.items():
            loading = is_loading_module
            is_loading_module = False
            terminal =  'Thioesterase' in module_domains.keys()
            domains_present = module_domains.keys()
            has_acp = 'ACP' in domains_present or 'PCP' in domains_present
            has_at = 'AT' in domains_present or 'CAL' in domains_present
            if not accept_in_trans_modules and not (has_at and has_acp): #
                print(f'\t\tModule with no {"AT" if has_acp else "ACP"} skipped: {",".join(domains_present)} \n')
                continue
            if not has_acp and accept_in_trans_modules:
                # No acp means it's either invalid or trans. Assume in trans LM for now
                module = pks.models.TransModule(subunit=db_subunit_entry, loading=loading, terminal=terminal)
            else:
                module = pks.models.Module(subunit=db_subunit_entry, loading=loading, terminal=terminal)
            if save_to_db:
                module.save()
                module.buildDomains(module_domains, cyclic=False) #no cycle information
            print(f'\t\t{",".join(module_domains)}')
            #TODO the antismash_to_database_functions file checks for an assertion erro. No idea if this is necessary





if clear_db:
    print('deleting old database...')
    [cluster.delete() for cluster in pks.models.Cluster.objects.all()]
    print('database cleared')


#filelist = [x for x in filelist if 'CP020044' in x]
for file in filelist:
    record = SeqIO.read(file, "genbank")
    #print(record)
    #continue
    if save_to_db and pks.models.Cluster.objects.filter(mibigAccession=record.annotations['accessions'][0]).exists():
        # We assume that mibig has better information than this pipeline will, so if we find
        if debug:
            print('skipped: found file in the database')
        continue
    if save_to_db and pks.models.Cluster.objects.filter(genbankAccession=record.annotations['accessions'][0]).exists():
        if debug:
            print('skipped: found file in the database')
        continue
    cluster = pks.models.Cluster(
        genbankAccession=record.annotations['accessions'][0],
        mibigAccession='', #does not have one because this is not mibig
        description='', # Problem: description is the known compound name, usually, but we don't have that here.
        sequence=record.seq,
        knownProductSmiles='',
        knownProductMCS=''
    )
    #print('\n\n\nprocessing {}'.format(file))
    #print(record.seq)
    #print(record)
    print(f'processed file {file}')
    subunits = read_subunit(record)
    if len(subunits) == 0:
        if debug:
            print('\tFile skipped due to lack of valid subunits.')
        continue
    if all(subunit['modules'] is None or len(subunit['modules']) == 0 for subunit in subunits):
        if debug:
            print('File skipped because all found subunits had no valid modules.')
        continue
    if save_to_db:
        cluster.save()
    process_subunits(subunits, cluster)
