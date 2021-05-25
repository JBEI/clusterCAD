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
import warnings

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


def db_save(item):
    if save_to_db:
        item.save()

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
            subunit_info['modules'] = antismash_funcs.processSubunitModules(feature.qualifiers['sec_met'], accept_in_trans_modules=accept_in_trans_modules)
            if debug:
                print('\t\t',feature.qualifiers['sec_met'])
                for order, module in subunit_info['modules'].items():
                    printable = ','.join(module.keys())
                    print(f'\t\t{order}: {printable}')
                print('-' * 200)
        elif pks_subtype == 'NRPS/PKS subtype: Type I Iterative PKS': 
            # Also process iterative PKS, with a note in the description
            subunit_info['description'] = 'Iterative PKS: ' + subunit_info['description']
            subunit_info['modules'] = antismash_funcs.processSubunitModules(feature.qualifiers['sec_met'], accept_in_trans_modules=accept_in_trans_modules)
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
            # No modules, could be standalone, but not supported for now.
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
        db_save(db_subunit_entry)
        # A list of bad modules that might be standalones.
        db_potential_standalones = []
        db_modules = []
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
            if debug:
                print('processing module {module_index}:', module_domains)
            if not accept_in_trans_modules and not (has_at and has_acp):
                print(f'\t\tModule with no {"AT" if has_acp else "ACP"} skipped: {",".join(domains_present)} \n')
                continue
            if not (has_at or has_acp): # invalid even for in trans LM
                print(f'\t\tModule with no AT or ACT skipped: {",".join(domains_present)} \n')
                continue
            if not has_acp and accept_in_trans_modules:
                # No acp means it's either invalid or trans. Assume in trans LM for now
                standalone = pks.models.DomainContainingStandalone(
                    cluster=cluster,
                    name=db_subunit_entry.name,
                    start=db_subunit_entry.start,
                    stop=db_subunit_entry.stop,
                    sequence=db_subunit_entry.sequence,
                    genbankAccession=db_subunit_entry.genbankAccession,
                )
                """Note: A standalone is at the level of a gene/subunit even though the
                domain containing aspect represents that of a module. Since we
                don't know whether there will be valid modules unil we process every single
                entry from "modules", then we will need to create these standalones in the
                meantime so that, if there end up being no valid modules, we save the
                standalones instead."""
                db_potential_standalones.append(standalone)
                # Save the created object as soon as it's made since it doesn't behave
                # like a python object and will cause bugs later if not.
                standalone.save()
            else:
                module = pks.models.Module(subunit=db_subunit_entry, loading=loading, terminal=terminal)
                # Save the created object as soon as it's made since it doesn't behave
                # like a python object and will cause bugs later if not.
                module.save()
                db_modules.append(module)
                #print(f'\t\t{",".join(module_domains)}')
        if save_to_db:
            """
            Save to the database our results. Since you can't have a standalone gene
            and a PKS Module in the same gene, we can only accept one. PKS Modules take
            precedence here and possible detected standalone genes will be saved otherwise.
            """
            # If valid modules exist, that takes precedence
            if len(db_modules) != 0:
                print('saving modules:', db_modules)
                for standalone in db_potential_standalones:
                    # save modules to db if there are any valid ones
                    standalone.delete()
                for module in db_modules:
                    module.buildDomains(module_domains, cyclic=False) #no cyclization information
            elif len(db_potential_standalones) != 0:
                print('saving standalones:', db_potential_standalones)
                if len(db_potential_standalones) != 1:
                    warnings.warn('Multiple standalones genes possible within the same gene.')
                # Delete subunit, since our subunit info is now in the standalone.
                db_subunit_entry.delete()
            else:
                print('subunit no good. Bye.')
                db_subunit_entry.delete()
            #TODO the antismash_to_database_functions file checks for an assertion erro. No idea if this is necessary

def reorder_subunits(cluster):
    """
    Attempts to predict the relative ordering of subunits inside
    a cluster using information about loading modules and thioesterases.
    Predictions of loading modules are dependent on in-trans or
    non-KS-AT-ACP modules, while the ends are predicted by thioesterases.
    """

    #### HELPER FUNCTIONS ####
    def must_be_first(subunit):
        """
        Given a subunit, looks at its first module to see if it is
        necessarily the first subunit
        """
        if (len(subunit.modules()) == 0):
            return False
        first_module = subunit.modules()[0]
        first_domain = first_module.domains()[0]
        if isinstance(first_domain, pks.models.AT) or isinstance(first_domain, pks.models.ACP):
            return True
        return False


    def must_be_last(subunit):
        """
        Given a subunit, looks at its last module to see if it is
        necessarily the last subunit
        """
        if (len(subunit.modules()) == 0):
            return False
        subunit_modules = subunit.modules()
        last_module = subunit_modules[len(subunit_modules) - 1]
        last_module_domains = last_module.domains()
        last_domain = last_module_domains[len(last_module_domains) - 1]
        if isinstance(last_domain, pks.models.TE):
            return True
        return False

    def panic(subunits):
        """
        Panic and set the order of everything to -1
        because we're not sure.
        """
        for subunit in subunits:
            subunit.order = -1
            db_save(subunit)

    def reorder_modules(cluster):
        if all(subunit.order != -1 for subunit in cluster.subunits()):
            # order things
            ordered_subunits = sorted(cluster.subunits(), key=lambda sub: sub.order)
            cur_order = 0
            for subunit in ordered_subunits:
                for module in subunit.modules():
                    print(cur_order)
                    module.order = cur_order
                    cur_order += 1
                    db_save(module)

        else:
            # TODO possibly partial ordering? but not the goal right now
            print('nothing going on here')
            pass
    #### END HELPER FUNCTIONS ####


    subunits = cluster.subunits()
    if len(subunits) not in [1, 2, 3]:
        return
    if len(subunits) == 1:
        subunits[0].order = 0
        subunits[0].save()
    if len(subunits) in [2, 3]:
        uncertain_subunits = set(subunits)
        first_subunit = None
        last_subunit = None
        for subunit in subunits:
            if must_be_first(subunit):
                if first_subunit is not None:
                    # Two possible "first" subunits? Something illegal.
                    return panic(subunits)
                first_subunit = subunit
                uncertain_subunits.remove(subunit)
            if must_be_last(subunit):
                lm_is_end_module = (first_subunit is subunit) and len(subunit.modules()) < 2
                if last_subunit is not None or lm_is_end_module:
                    # can't be both LM and last subunit
                    return panic(subunits)
                last_subunit = subunit
                if subunit in uncertain_subunits:
                    uncertain_subunits.remove(subunit)
        if len(uncertain_subunits) == 1:
            remaining = uncertain_subunits.pop()
            remaining.order = 0 if first_subunit is None else 1
            db_save(remaining)
        else:
            # Don't know what the order of the rest of the subunits are.
            panic(uncertain_subunits)
        if last_subunit:
            last_subunit.order = len(subunits) - 1
            db_save(last_subunit)
        if first_subunit:
            first_subunit.order = 0
            db_save(first_subunit)

    if all(subunit.order != -1 for subunit in cluster.subunits()):
        ordered_subunits = [sub.name for sub in sorted(cluster.subunits(), key=lambda sub: sub.order)]
        # Handles module ordering and loading modules
        cluster.reorderSubunits(ordered_subunits)
        print([(subunit.modules(), subunit.order) for subunit in subunits])




def filter_cluster(cluster):
    """
    Filters cluster if, throughout ALL of the mudoles in the subunit, we don't
    have a _single_ KS, AT, or ACP.
    """
    domains_present = set()
    for subunit in cluster.subunits():
        for module in subunit.modules():
            domains_present.update([repr(domain) for domain in module.domains()])
        has_ks = 'KS' in domains_present
        has_at = 'AT' in domains_present or 'CAL' in domains_present
        has_acp = 'ACP' in domains_present or 'PCP' in domains_present
        if has_ks and has_at and has_acp:
            return
    print(f'bad cluster deleted: {cluster.genbankAccession, cluster.description}')
    print('contained', domains_present)
    cluster.delete()




if clear_db:
    print('deleting old database...')
    [cluster.delete() for cluster in pks.models.Cluster.objects.all()]
    print('database cleared')


#filelist = [x for x in filelist if '01854385' in x]
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
    db_save(cluster)
    process_subunits(subunits, cluster)
    if save_to_db:
        # Only clear data if we saved it to begin with
        filter_cluster(cluster)
    reorder_subunits(cluster)
