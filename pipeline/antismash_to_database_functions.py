#!/usr/bin/python

import os,sys
import glob
import json
import re
from collections import OrderedDict

from Bio import SeqIO

allowed_domains = ['KS', 'AT', 'KR', 'DH', 'ER', 'ACP', 'Thioesterase', 
                      'cMT', 'oMT', 'CAL', 'PCP', 
                      'Heterocyclization', 'AMP-binding', 
                      'Condensation_Starter',
                      'Condensation_DCL', 'Condensation_LCL',
                      'PKS_Docking_Nterm', 'PKS_Docking_Cterm']
sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

def mibigSubtypes(filepath):
    '''Takes as input a file path and outputs set of all PKS subtypes that 
       appear in all MIBiG JSON files that are found along that path. 
    '''
    filelist = glob.glob(os.path.join(filepath, '*.json'))
    subtypes = []
    for filename in filelist:
        with open(filename) as jsonfile:
            jsondata = json.load(jsonfile)
        try:
            subtypes.extend(jsondata['general_params']['Polyketide']['pks_subclass'])
        except KeyError:
            pass
    return set(subtypes)

def filterModularTypeI(filepath, validset):
    '''Takes as input a file path and a set of str indicating what
       MIBiG PKS subunit annotations should be considered as valid type I
       modular PKSs. Outputs list of MIBiG accession numbers that have an
       annotated PKS subtype that is consistent with the valid set.
    '''
    assert isinstance(validset, set)
    filelist = glob.glob(os.path.join(filepath, '*.json'))
    accessions = []
    for filename in filelist:
        with open(filename) as json_file:
            json_data = json.load(json_file)
        try:
            if len(validset.intersection( \
              set(json_data['general_params']['Polyketide']['pks_subclass']))) > 0:
                accession = os.path.basename(filename).strip('.json')
                accessions.append(accession)
        except KeyError:
            pass
    return accessions

def processSubunitModules(sec_met): 
    '''Takes as input annotated features for a PKS subunit following 
       antiSMASH analysis and returns a dict of OrderedDict objects 
       corresponding to the annotated modules contained by the subunit. The 
       keys of the subunit OrderedDict are the index of the module within the
       subunit. This function assumes that feature.type=='CDS' and that 
       feature.qualifiers has the key 'sec_met'. If the first entry of 
       'sec_met' is not 'Type: t1pks' then nothing is returned.
    '''
    # Initialize dict for representation of the subunit
    # keys: module number
    # values: OrderedDict of domains in module. Within the OrderedDict, 
    #         the key is the domain name, while the value is a length one or 
    #         two list where the first element is a dictionary 
    #         {start: value, stop: value} and the optional second element is 
    #         a specificity dict if applicable
    subunit = {}
    
    module_index = 0 # track current module
    module_domains = [] # list of domains in module
    old_module_domains = [] # pre-initialize in case module starts with a domain
                            # that is expected to end the module
    
    # This is how domains appear in sec_met:
    # ['PKS_AT', 'PKS_KS', 'PKS_KR', 'PKS_DH', 'PKS_ER', 'ACP', 'Thioesterase']
    # Iterate over the entries in sec_met, and add them to the module_domains list 
    for entry in sec_met:    
        # Split entry into a list
        entrysplit = [item.strip() for item in entry.split(';') if item != '']
        # Split part of entry that is expected to describe catalytic domain
        domainsplit = entrysplit[0].split()
        if not (' '.join(domainsplit[:2]) == 'NRPS/PKS Domain:' and len(domainsplit) > 2):
            # Should be a predicted PKS Domain
            continue
        special_cases = {
            'PKS_Docking_Nterm': 'PKS_Docking_Nterm',
            'PKS_Docking_Cterm': 'PKS_Docking_Cterm',
            'CAL_domain': 'CAL',
            # Assume 'DH2' and 'DHt' are the same as 'DH' 
            'PKS_DH2': 'DH',
            'PKS_DHt': 'DH'
        }
        domaintype = domainsplit[2]
        if domaintype in special_cases:
            domaintype = special_cases[domaintype]
        elif domaintype[:4] == 'PKS_':
            # If not a special case and there is a leading PKS, trim it
            domaintype = domaintype[4:]
        # These are the catalytic domains that ClusterCAD wil recognize
        if domaintype not in allowed_domains:
            print('\tIgnoring domain type: %s' %(domaintype))
            # Break out of for loop and stop looking for additional catalytic domains if 
            # we encounter a domain that we don't recognize
            # We end up excluding any subunit that has a non-recognized catalytic domain
            break    
        # Get the boundaries of the catalytic domain
        boundaries = [int(bound) for bound in re.sub(r'[().]', '', domainsplit[3]).split('-')]

        # AntiSMASH seems to count from 0 for start positions
        # but 1 for stop positions, as in BioPython
        # so we add 1 to the start start here
        boundaries[0] += 1
        
        #print(domaintype)
        # Here, we add each domain to a list, which will be converted to an OrderedDict
        # Include substrate and stereospecificity annotations for CAL, AT, and KR domains respectively
        if domaintype in ['AT', 'KR', 'CAL']:
            notesdict = {}
            for note in entrysplit[1:]:
                item = note.split(': ')
                notesdict[item[0]] = item[1]
            module_domains.append((domaintype, 
                                   [{'start': boundaries[0], 'stop': boundaries[1]}, notesdict]))
 
        # End of the module has been reached of the domain is 'ACP' or 'PCP
        elif domaintype in ['ACP', 'PCP']:
            module_domains.append((domaintype, 
                                   [{'start': boundaries[0], 'stop': boundaries[1]}]))
            domains_present = [d[0] for d in module_domains]
            # Make sure every module has an AT or CAL, or else it isn't valid and should be ignored
            # This means it will be excluded from the subunit, which makes sense since we can't 
            # really perform a polyketide chain extension without an AT
            if 'AT' in domains_present or 'CAL' in domains_present:            
                subunit[module_index] = OrderedDict(module_domains)
                old_module_domains = module_domains
                module_index += 1
            else:
                # UPDATE: keeping ACP/PCPs with no AT and CAL, in the case 
                # of in trans loading modules.
                print(f'added in trans loading module {",".join(module[0] for module in module_domains)}')
                subunit[module_index] = OrderedDict(module_domains)
                old_module_domains = module_domains
                module_index += 1

                old_module_domains = []
            module_domains = []
        # These domains may come after the ACP or PCP, so if they are encountered, we add
        # them to previous module and keep going forward
        elif domaintype in ['Thioesterase', 'PKS_Docking_Cterm', 'Condensation_LCL']:
            # Overwrite previous subunit, or else will have duplicate entries
            old_module_domains.append((domaintype, 
                                       [{'start': boundaries[0], 'stop': boundaries[1]}]))
            subunit[module_index - 1] = OrderedDict(old_module_domains)
            module_domains = []
        #for all other domains, assume the same
        elif domaintype in allowed_domains:
            module_domains.append((domaintype, 
                                   [{'start': boundaries[0], 'stop': boundaries[1]}]))
    # Don't throw away any leftovers if they contain AT, KS, or CAL, for 
    # in trans loading modules.
    if module_domains and any(domain[0] in ['AT', 'KS', 'CAL'] for domain in module_domains):
        subunit[module_index] = OrderedDict(module_domains)
        print(f'added in trans loading module {",".join(module[0] for module in module_domains)}')
    return subunit

def processClusterSeqRecord(record):
    '''Takes as input a SeqRecord read in from an annotated MIBiG .embl file 
        and outputs list of sublists containing PKS data from that record. 
        Output is a list where each entry is a sublist that represents 
        a subunit or standalone.
    '''
    # Get list to hold information about all genes in the record
    gene_data = []

    
    for feature in record.features:
        # These are the features we are interested in
        if feature.type == 'CDS' and 'protein_id' in feature.qualifiers.keys() and \
                ('gene' in feature.qualifiers.keys() or 'product' in feature.qualifiers.keys()): 
            # This gets the location of the feature
            location = feature.location
            # General information about gene
            if 'product' in feature.qualifiers.keys():
                description = feature.qualifiers['product'][0]
            if 'gene' in feature.qualifiers.keys(): 
                genename = feature.qualifiers['gene'][0]
            elif 'product' in feature.qualifiers.keys():
                genename = feature.qualifiers['product'][0]
            gene_data.append([feature.qualifiers['protein_id'][0],
                              genename,
                             ])
            # Feature may not be a PKS module and therefore may not have have subunits 
            # (this will be overwritten if it does have subunits)
            subunit_modules = None
            # Information if gene is PKS subunit
            if 'sec_met' in feature.qualifiers.keys() and len(feature.qualifiers['sec_met']) > 3:
                if feature.qualifiers['sec_met'][3] in \
                  ['NRPS/PKS subtype: Type I Modular PKS', 
                   'NRPS/PKS subtype: PKS-like protein',
                   'NRPS/PKS subtype: PKS/NRPS-like protein',
                   'NRPS/PKS subtype: Hybrid PKS-NRPS']:
                    subunit_modules = processSubunitModules(feature.qualifiers['sec_met'])

            # Append description and position of gene within nucleotide sequence
            gene_data[-1].extend([description, [location.start.position + 1, location.end.position]])

            # Subunit information (if no subunit information, assumed to be a standalone enzyme)
            if subunit_modules:
                gene_data[-1].append(subunit_modules)

            # Amino acid sequence
            gene_data[-1].append(feature.qualifiers['translation'][0])

    return gene_data

def checkModuleValidity(modulelist):
    '''Function that makes sure module specified in MIBiG JSON file is valid,
       that is to say, make sure that it contains KS, AT, and ACP or PCP.
       These domains have the following possible annotations:
       ['KS', 'AT', 'T']
       ['Ketosynthase', 'Acyltransferase', 'Thiolation (ACP/PCP)']
    '''
    atcheck = len(set(['AT', 'CAL', 'Acyltransferase']).intersection(set(modulelist)))
    acpcheck = len(set(['ACP', 'PCP', 'T', 'Thiolation (ACP/PCP)']).intersection( \
      set(modulelist)))
    if atcheck and acpcheck:
        return True
    else:
        return False

def enterCluster(cluster, clusterrecord, mibigfile):
    '''Takes in a reference to a Django cluster object, a record opened from a
       antiSMASH embl file and corresponding MiBiG JSON file, then uses these 
       enters cluster information into Django database.
    '''

    # Get information about the gene
    gene_data = processClusterSeqRecord(clusterrecord)

    if len(gene_data) == 0:
        return

    # Initalize lists for subunits and standalones
    # We make two dictionaries because sometimes the subunit name in the MiBiG JSON files
    # is the gene name, e.g. eryA1, and sometimes it is the accession number, e.g. A0000000
    unordered_subunits = {}
    unordered_subunits_alt = {}
    standalones = []
     
    # Recall that each entry in gene_data is a list
    # [protein id, gene, product, [location start, location end], 
    #  subunit dict (optional), translation]
    
    #####################
    # Basic information #
    #####################
    
    counter = 1
    for gene in gene_data:
        geneid = gene[0].strip()
        genename = gene[1].strip()
        genedesc = gene[2].strip()
        genestart = gene[3][0]
        genestop = gene[3][1]
        genetranslation = gene[-1].strip()

        # Just use length of gene_data to differentiate between standalones and subunits
        if len(gene) == 6:
            # We do this to take care of duplicated gene names
            if genename in unordered_subunits_alt.keys():
                genename = genename + '_' + str(counter)
                counter += 1
 
            # Get subunit data from gene
            genesubunitdata = gene[-2]
            # Here we use the two dictionary options to save the unordered subunits
            # Sometimes MIBiG uses geneid and sometimes it uses genename to reference subunits
            unordered_subunits[geneid] = (genename, genedesc, genestart, genestop,
                                            genesubunitdata, genetranslation)
            unordered_subunits_alt[genename] = (geneid, genedesc, genestart, genestop,
                                                genesubunitdata, genetranslation)
        else:
            # Standalones lack subunit and orphan entries
            assert len(gene) == 5, gene
            
    #####################
    # CREATE STANDALONE #
    #####################

#            pks.models.Standalone(cluster=cluster)
#            standalones.append(pks.Standalone(geneid, genename, genedesc, 
#                                              genestart, genestop, genetranslation))
    
    #########################################
    # JSON file has cyclization information #
    #########################################

    # Get ordered version of subunits from corresponding JSON file
    with open(mibigfile) as json_file:
        mibig_data = json.load(json_file)
    
    # Get PKS cyclization information
    # this will be either 'Cyclic' or 'Linear'
    try:
        lin_cycl_pk = mibig_data['general_params']['Polyketide']['lin_cycl_pk']
        if lin_cycl_pk == 'Cyclic':
            cyclize = True
        elif lin_cycl_pk == 'Linear':
            cyclize = False
        else:
            raise Exception("lin_cycl_pk expected to be 'Cyclic' or 'Linear'.")
    except KeyError:
        cyclize = False
            
    #############################################
    # Case 1: Use JSON file subunit information #
    #############################################
        
    # Note that all gene data has now been processed, want to reprocess to get right ordering 
    # We strip out subunits that have invalid modules
    try:
        ordered_subunits = []
        for subunit in mibig_data['general_params']['Polyketide']['mod_pks_genes']:
            subunit_name = re.sub(r'\s+', '', subunit['mod_pks_gene'])
            subunit_modules = subunit['pks_module']

            valid_subunit = True
            # This checks if the module is valid
            for module in subunit_modules:
                if not check_json_module_validity(module['pks_domains']):
                    valid_subunit = False
            if valid_subunit:
                ordered_subunits.extend(subunit_name.split(','))
            else:
                # Loop is broken once first invalid subunit is encountered
                break
        # If no valid subunits, then raise exception to use alphabetical ordering
        if len(ordered_subunits) == 0:
            raise Exception

        # This makes sure the subunit accession naming is consistent
        # The purpose of these two 'if' statements is because there may be cases 
        # in the MiBiG JSON file where the name of the gene is for example, 
        # 'eryA1, A000000' and we want to keep consistant naming
        if len(ordered_subunits[0]) >= 8:
            ordered_subunits = [entry for entry in ordered_subunits if len(entry) >= 8]
        if len(ordered_subunits) > 1:
            if len(ordered_subunits[1]) >= 8:
                ordered_subunits = [entry for entry in ordered_subunits if len(entry) >= 8]
        # This is because the accession number under which the gene is recorded sometimes
        # has a version number, and sometimes does not
        if len(ordered_subunits[0].split('.')) == 1 and len(ordered_subunits[0]) == 8:
            ordered_subunits = [entry + '.1' for entry in ordered_subunits]

        # Check if subunit is in either dictionary
        for isubunit,subunit in enumerate(ordered_subunits):
            if subunit not in set(list(unordered_subunits.keys()) + \
              list(unordered_subunits_alt.keys())):
                print('Missing subunit: "%s"' %(subunit))
                for gene in mibig_data['general_params']['Polyketide']['mod_pks_genes']:
                    if gene['mod_pks_gene'] == subunit:
                        module = gene['pks_module']
                        for entry in module:
                            print(entry['pks_domains'])
                print(unordered_subunits.keys())
                return

        # Determine whether to use standard or alternative dict
        if len(ordered_subunits[0]) >= 8:
            alt = False
        else:
            alt = True
    
    #########################################
    # Case 2: Alphabetical subunit ordering #
    #########################################
    
    # Just use unordered gene order if the gene ordering is not already in the JSON file 
    except Exception:
        ordered_subunits = list(unordered_subunits_alt.keys())
        ordered_subunits.sort()
        alt = True

    ####################################
    # This does the subunit reordering #
    ####################################
    # This is just to make sure loading module is assigned correctly
    modules_seen = 0
    for subunit_key in ordered_subunits:
        if not alt:
            subunitdata = unordered_subunits[subunit_key]
        else:
            subunitdata = unordered_subunits_alt[subunit_key]
     
        if not alt:
            subunit = pks.models.Subunit(cluster=cluster,
                                         genbankAccession=subunit_key,
                                         name=subunitdata[0],
                                         start=subunitdata[2],
                                         stop=subunitdata[3],
                                         sequence=subunitdata[-1])
            subunit.save()
        else:
            subunit = pks.models.Subunit(cluster=cluster,
                                         genbankAccession=subunitdata[0],
                                         name=subunit_key,
                                         start=subunitdata[2],
                                         stop=subunitdata[3],
                                         sequence=subunitdata[-1])
            subunit.save()
        # These are the modules for the subunit
        moduledata = subunitdata[-2]
        # We lump in the loading didomain and TE on the first and last modules respectively
        modulekeys = list(moduledata.keys())
        #index of the module
        imodule = 0
        while imodule < len(modulekeys):
            # Get info
            keys = list(moduledata[modulekeys[imodule]].keys())
            values = moduledata[modulekeys[imodule]].values()
            # Process info according to loading or not
            if modules_seen == 0:
                loading = True                
            else: 
                loading = False
            moduledict = OrderedDict([(k,v) for k,v in zip(keys,values)])
            # Determine whether module is terminal or not
            if 'Thioesterase' in list(moduledata[modulekeys[imodule]].keys()):
                terminal = True
            else:
                terminal = False
            imodule += 1
            modules_seen += 1
            try:
                # This is to make sure we don't add subunits with invalid modules
                # The check for errors here is to compare against the known product
                domains_present = moduledict.keys()
                if 'ACP' in domains_present or 'PCP' in domains_present:
                    if 'AT' in domains_present or 'CAL' in domains_present:
                        module = pks.models.Module(subunit=subunit, 
                                                   loading=loading, terminal=terminal)
                        module.save()
                        module.buildDomains(moduledict, cyclic=cyclize)
                        continue
                print(f'invalid module with {domains_present} not saved')
            except AssertionError as e:
                print(moduledict)
                print(type(e).__name__, e.args, subunit + ' ' + subunitdata[1])
                raise Exception(type(e).__name__, e.args, subunit + ' ' + subunitdata[1])
                break
