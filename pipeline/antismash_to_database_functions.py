#!/usr/bin/python

import os,sys
import glob
import json
from collections import OrderedDict

from Bio import SeqIO

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

def buildDomain(module, domain_data):
    '''Takes a given value from a dictionary and builds the domain'''
    
    default_type = domain_data[0]
    domain_subtype = domain_data[1]
    start = domain_data[2]
    stop = domain_data[3]
    specificity = domain_data[4]
        
    # Process antiSMASH PKS and NRPS domain type into more commonly known type abbreviations
    
    # Trim off leading "PKS_" for some PKS domains
    if default_type.split("_")[0] == "PKS":
        if str(default_type) in ['PKS_Docking_Nterm', 'PKS_Docking_Cterm']:
            domain_type = default_type
        elif default_type.split("_")[1] == "PP":
            domain_type = 'ACP'
        else:
            # Assume 'DH2' and 'DHt' are the same as 'DH' 
            domain_type = default_type.split('_')[-1].replace('DHt', 'DH').replace('DH2', 'DH')
    
    # Methyltransferases have subtypes: oMT, cMT, nMT. We want the subtype to be the type.
    elif default_type == 'MT':
        domain_type = domain_subtype
    
    # Condensation domains have subtypes: Starter, LCL, DCL, Dual, Cglyc, etc. We want the subtype to be the type.
    # Some however don't have subtypes, so just use Condensation if subtype is not listed.
    elif default_type == 'Condensation': 
        if domain_subtype != "N/A":
            domain_type = domain_subtype
        else:
            domain_type = default_type
    
    elif default_type == 'CAL_domain':
        domain_type = 'CAL'
    elif default_type == 'AMP-binding':
        domain_type = 'A'
    elif default_type == 'Epimerization domain':
        domain_type = 'E'
    elif default_type == 'Epimerization':
        domain_type = 'E'
    elif default_type == 'Heterocyclization':
        domain_type = 'Cy'
    elif default_type == 'TD':
        domain_type = 'R'
    elif default_type == 'A-OX': # substrate is glycine but effectively chain terminating
        domain_type = 'AOX'
    elif default_type in ['NRPS-COM_Nterm', 'NRPS-COM_Cterm']:
            domain_type = default_type
    elif default_type == 'PP-binding':
        domain_type = 'ACP'
    else:
        domain_type = default_type
    
    # Make sure the domain of interest is recognized (i.e. a model exists for it)
    # If domain is not recognized, no domain object is generated.  
    
    if domain_type in ['KS', 'AT', 'KR', 'DH', 'ER', 'ACP', 'Thioesterase', 
                              'cMT', 'oMT', 'CAL', 'PCP', 
                              'Cy', 'A', 'Condensation', 'Condensation_Starter',
                              'Condensation_DCL', 'Condensation_LCL', 'Condensation_Dual', 'Cglyc',
                              'PKS_Docking_Nterm', 'PKS_Docking_Cterm',
                              'AOX', 'E', 'F', 'nMT', 'R', 'X',
                              'NRPS-COM_Nterm', 'NRPS-COM_Cterm']:
        #print("Domain recognized, building domain. ")
        
        # Directly build domains that have no special conditions besides active/inactive
        if domain_type in ['KS', 'DH', 'ER', 'cMT', 'oMT', 'PCP', 'Cy', 'E', 'nMT', 'F', 'AOX', 'X']:
            if domain_type in ['DH', 'ER', 'cMT', 'oMT', 'nMT', 'F', 'AOX', 'X']:
                # getattr(pks.models, domain_type) will fetch domain_type from pks.models
                # and replace itself with the target domain class object
                newDomain = getattr(pks.models, domain_type)(module=module, start=start, stop=stop, active=True)
            else:
                newDomain = getattr(pks.models, domain_type)(module=module, start=start, stop=stop)
            newDomain.save()
        
        # Build domains that have special conditions like substrate specificity or stereoselectivity
        if domain_type == 'ACP':
            # No default ACP substrate, but adjusted to mal (manually adjustable) if no AT detected (N-ACP) 
            substrate = 'None'
            if len(pks.models.AT.objects.filter(module=module)) == 0:
                substrate = 'mal'
            newDomain = pks.models.ACP(module=module, start=start, stop=stop, substrate=substrate)
            newDomain.save()
        if domain_type in ['CAL', 'AT']:
            substrate = specificity[0].split()[1]
            if (substrate not in pks.models.starters) and (substrate not in pks.models.extenders):
                #print("Substrate not recognized, assigning as malonyl. ")
                substrate = 'mal'
            newDomain = getattr(pks.models, domain_type)(module=module, start=start, stop=stop, substrate=substrate)
            newDomain.save()
        if domain_type == 'KR':
            activity = specificity[0].split()[1]
            if activity == 'active':
                active = True
            else:
                active = False
            KR_type = specificity[1].split()[2]
            if KR_type == '(unknown)':
                KR_type = 'U'
            newDomain = pks.models.KR(module=module, start=start, stop=stop, active=active, type=KR_type)
            newDomain.save()
        if domain_type == 'A':
            substrate = specificity[0].split()[1]
            # if substrate is invalid or not implemented, assign as alanine (glycine causes issues with E domains)
            if (substrate not in pks.models.starters) and (substrate not in pks.models.extenders):
                #print("Substrate not recognized, assigning as ala. ")
                substrate = 'ala'
            # A redox type defaults to none; can be manually adjusted to have oxidation/reduction loop out in Cy's presence
            A_type = 'none'
            newDomain = pks.models.A(module=module, start=start, stop=stop, substrate=substrate, type=A_type)
            newDomain.save()
        if domain_type in ['Condensation', 'Condensation_Starter', 'Condensation_LCL', 'Condensation_DCL', 
                           'Condensation_Dual', 'Cglyc']:
            if domain_type == 'Condensation':
                C_type = 'Unspecified'
            elif domain_type == 'Cglyc':
                C_type = 'Glycopeptide'
            else:
                C_type = domain_type.split('_')[1]
            newDomain = pks.models.C(module=module, start=start, stop=stop, type=C_type)
            newDomain.save()
        if domain_type == 'R':
            R_type = 'alcohol'
            newDomain = pks.models.R(module=module, start=start, stop=stop, active=True, type=R_type)
            newDomain.save()
        if domain_type == 'Thioesterase':
            # Cyclization is FALSE until corrections files kick in
            newDomain = pks.models.TE(module=module, start=start, stop=stop, cyclic=False, ring=0)
            newDomain.save()
        module.setLoading()
    
def checkModuleValidityAndReductiveDomainActivity(module):
    '''Examines a generated module for validity and deletes it if the module is invalid. 
       Also checks reductive domains for potential dependency issues affecting activity.
    '''
    
    domain_list = []
    for newdomain in module.domains():
        domain_list.append(newdomain.__repr__())

    #print(domain_list)
    
    # Each list contains elements that must be present in a module (excluding optional domains) for it to be a valid module
    # all() checks if all elements of specified list is in domain_list in some order
    # AOX also acts as A domain (A-OX)
    # Lone ACP allowed for N-ACPs
    # KS is checked for before an operation is applied for non-loading PKS modules
    if all(i in domain_list for i in ['AT','ACP']) or \
        all(i in domain_list for i in ['AT','PCP']) or \
        all(i in domain_list for i in ['A','ACP']) or \
        all(i in domain_list for i in ['A','PCP']) or \
        all(i in domain_list for i in ['CAL']) or \
        all(i in domain_list for i in ['C','A','PCP']) or \
        all(i in domain_list for i in ['Cy','A','PCP']) or \
        all(i in domain_list for i in ['C','AOX','PCP']) or \
        domain_list == ['ACP']:        
        
        # Now check reductive domains for potential dependency issues affecting activity
        # If there is a DH and no KR, DH should be inactivated
        # If there is a ER and no KR (or no DH), ER should be inactivated
        if ('DH' in domain_list) and ('KR' not in domain_list):
            print("Warning: DH present without KR. Inactivating DH. ")
            dh = pks.models.DH.objects.filter(module=module)[0]
            print(dh)
            dh.active = False
            dh.save()
        if ('ER' in domain_list) and ('DH' not in domain_list):
            print("Warning: ER present without DH. Inactivating DH. ")
            er = pks.models.ER.objects.filter(module=module)[0]
            er.active = False
            er.save()
        if ('ER' in domain_list) and ('KR' not in domain_list):
            print("Warning: ER present without KR. Inactivating DH. ")
            er = pks.models.ER.objects.filter(module=module)[0]
            er.active = False
            er.save()
            
        return True
    else:
        module.delete()
        print("***Incomplete module detected: Module deleted. ***")
        return False

def reverseModuleOrder(newSubunit):
    '''Generates and reverses a list of module indices within the subunit 
    and reassigns module order accordingly'''
    
    old_module_order = []
    for newmodule in newSubunit.modules():
        old_module_order.append(newmodule.order)
    new_module_order = list(reversed(old_module_order))
    #print(old_module_order)
    #print(new_module_order)

    index = 0
    for newmodule in newSubunit.modules():
        order = new_module_order[index]
        newmodule.order = order
        newmodule.save()
        index += 1

def reverseSubunitOrder(cluster):
    '''Generates and reverses a list of subunit names and reorders subunits accordingly'''
    
    old_subunit_order = []
    for newsubunit in cluster.subunits():
        old_subunit_order.append(newsubunit.name)
    new_subunit_order = list(reversed(old_subunit_order))
    cluster.reorderSubunits(new_subunit_order)
    
        
def buildCluster(cluster, record):
    '''Aggregates subunit, module, and domain information for a cluster 
    from an antiSMASH output file and builds the cluster'''
    
    #################### STAGE 1: PARSE ALL DATA FROM ANTISMASH FILE ######################

    # First compile dictionaries of subunits, modules, and domains from antiSMASH generated features
    # All relevant data will be stored in dictionaries for later access
    
    subunits_dict = {} 
    # key/value format: {gene/locus_tag : 0-gbk_acc, 1-name, 2-start, 3-stop, 4-seq, 5-sense, 6-alternative name(gene)}
    # key (either gene/locus tag) is what all modules and domains are linked to
    # name is same as key; alternative name is gene (common name) if available, otherwise N/A
    
    # counter used to append numbers if identical subunit common names (alt name) are detected
    same_name_counter = 1
    
    modules_dict = {} 
    # key/value format: {locus_tags : [(0-loading, 1-terminal, 2-list of domains), (load2, term2, domainlist2)]}
    # data for each individual module is stored in a list because all modules within a subunit share the same locus_tags
    
    domains_dict = {} 
    # key/value format: {domain_id : [0-default_type, 1-domain_subtype, 2-start, 3-stop, 4-specificity]}
    # specificity is a list with format [KR activity, KR type] for KR and [predict. substrate, predict. substrate,..] for AT/A
    
    # Iterate through all features, collecting and parsing data (qualifiers) required to construct objects
    for feature in record.features:
        
        # Clusters can use "locus_tag" or "gene" to identify subunits so we need to check both. 
        # Locate subunits using CDS and locus_tag/gene -- CDS only exist for genes/subunits
        if (feature.type == 'CDS') and \
            (("locus_tag" in feature.qualifiers) or ("gene" in feature.qualifiers) or ("product" in feature.qualifiers)):
            
            #print(feature.qualifiers["locus_tag"])
            
            subunit_gbk_accession = feature.qualifiers["protein_id"][0]
            
            # Subunit_alt_name saves "gene" to be used as the more common subunit name (e.g. eryAI) if available
            # Check for "gene" first so:
            # subunit_name can be overwritten by "locus_tag" (which module/domain could depend on) if applicable
            # AND gene overwrites N/A if available
            subunit_alt_name = "N/A"
            if "product" in feature.qualifiers:
                subunit_name = feature.qualifiers["protein_id"][0]
                subunit_alt_name = feature.qualifiers["product"][0]
            if "gene" in feature.qualifiers:
                subunit_name = feature.qualifiers["gene"][0]
                subunit_alt_name = feature.qualifiers["gene"][0]
            if "locus_tag" in feature.qualifiers:
                subunit_name = feature.qualifiers["locus_tag"][0]
            

            # start + 1 compensates for -1 property (seq->python) of the Biopython parser for sequence indices.
            subunit_start = feature.location.start.position + 1
            subunit_stop = feature.location.end.position
            subunit_sequence = feature.qualifiers["translation"][0]

            # Check if this subunit is on the complement (antisense) based on location [start,stop](+/-). 
            # If so, we need to reverse module order later because of reversed parsing order
            subunit_sense = str(feature.location).split("]")[1]
            
            # Check if this subunit shares a common name (alt_name) with another subunit in the same cluster
            # If so, append a "_1", "_2", etc. to distinguish them
            for s in subunits_dict.values():
                if s[6] == subunit_alt_name:
                    subunit_alt_name = "_".join([subunit_alt_name, str(same_name_counter)])
                    same_name_counter += 1
            
            subunits_dict[subunit_name] = (subunit_gbk_accession, subunit_name, subunit_start, 
                                           subunit_stop, subunit_sequence, subunit_sense, subunit_alt_name)
            
        # Locate modules using aSModule
        elif (feature.type == 'aSModule'):
            
            # Ensure this module only belongs to one subunit
            assert len(feature.qualifiers["locus_tags"]) == 1, "This module belongs to multiple subunits..."
            module_subunit = feature.qualifiers["locus_tags"][0]

            # If the qualifier starter_module is in the module's keys, then it is a starter/loading module
            if "starter_module" in feature.qualifiers: 
                module_loading = True
            else:
                module_loading = False

            # Modules with a TE or R or AOX are definitely terminating. 
            # Modules are reordered later so 
            module_terminal = False

            # Add list of domains within module into module_domains list
            module_domains = []
            module_domains = feature.qualifiers["domains"]
                        
            # Since multiple modules can be in a subunit, store module data in a list
            if module_subunit not in modules_dict:
                mod_data = [(module_loading, module_terminal, module_domains)]
                modules_dict[module_subunit] = mod_data
            else:
                # appends tuple to list of module under subunit key
                modules_dict[module_subunit].append((module_loading, module_terminal, module_domains))
        
        # Locate domains using aSDomain
        elif (feature.type == 'aSDomain'):
            
            domain_id = feature.qualifiers["domain_id"][0]
            
            # Fetch aSDomain type as default_type
            default_type = feature.qualifiers["aSDomain"][0]

            if "domain_subtype" in feature.qualifiers:
                domain_subtype = feature.qualifiers["domain_subtype"][0]
            else:
                domain_subtype = "N/A"
            
            domain_start = str(int(feature.qualifiers["protein_start"][0]) + 1)   #to account for 1-index of sequences
            domain_stop = feature.qualifiers["protein_end"][0]
            
            if "specificity" in feature.qualifiers:
                specificity = feature.qualifiers["specificity"] # specificity is a list
            else:
                specificity = "N/A"
            
            domains_dict[domain_id] = (default_type, domain_subtype, domain_start, domain_stop, specificity)
        
    #################### STAGE 2: CONSTRUCT SUBUNITS, MODULES, AND DOMAINS ######################
    
    #print("subunit keys")
    #print(subunits_dict.keys())
    #print("module keys")
    #print(modules_dict.keys())
    #print(modules_dict.values())
    
    # Now that we have identified all the subunits with modules in a cluster, create each subunit
    for subunit in subunits_dict.keys():
        
        subunitdata = subunits_dict.get(subunit)
        module_counter = 0
        
        #print(subunit)
        #print(type(subunitdata))
        #print("######")
        #print(subunit)
        
        #print(str(subunitdata[2])+"-"+str(subunitdata[3])) #location
        #print(subunitdata[5]) # sense
        
        # If available, use gene as subunit name since it's typically the common name
        #print(subunitdata[6])
        if subunitdata[6] != "N/A":
            subunit_name = subunitdata[6]
        else:
            subunit_name = subunitdata[1]
        
        #print(subunit_name)
        
        #print("Attempting to build Subunit %s." % subunit_name)
        newSubunit = pks.models.Subunit(cluster=cluster,
                             genbankAccession=subunitdata[0],
                             name=subunit_name,
                             start=subunitdata[2],
                             stop=subunitdata[3],
                             sequence=subunitdata[4], 
                             sense=subunitdata[5])
        newSubunit.save()
        #print("***Subunit %s built. ***" % subunit_name)

        # Now that each subunit has been generated, build each module
        # modulekey is the subunit name (locus_tags)
        # module is the list of modules+data within the subunit
        for modulekey in modules_dict.keys(): 
            
            module_list = modules_dict.get(modulekey)

            # Locate and build modules using locus_tags
            if modulekey == subunit:
                
                #print(modulekey)
                #print(len(module_list))
                subunit_module_counter = 0
                
                # Iterate through and build list of modules and domains (located at 2nd index of values) under a key
                while subunit_module_counter < len(module_list): 
              
                    #print("Attempting to build Module %d. " % module_counter)
                    newModule = pks.models.Module(subunit=newSubunit, 
                                               loading=module_list[subunit_module_counter][0], 
                                               terminal=module_list[subunit_module_counter][1])
                    newModule.save()
                    #print("***Module %d built. ***" % module_counter)

                    # Iterate through and build domains within this module
                    # Locate and build domains using domain_id as key
                    for domain in module_list[subunit_module_counter][2]:

                        #print(domain)
                        domain_data = domains_dict.get(domain)

                        #print("Attempting to build domain %s. " %(domain))
                        buildDomain(newModule, domain_data)
                        #print("Domain %s built. " %(domain))

                    #print("***All domains loaded into Module %d. ***" % module_counter)

                    # Check if the generated Module is acceptable by retrieving and checking newly generated domains
                    # If module is good, function returns True, if not then False. 
                    # Also checks for reductive domain dependency issues and adjusts activity accordingly
                    if checkModuleValidityAndReductiveDomainActivity(newModule):
                        module_counter += 1
                    else: 
                        module_counter += 0
                        
                    subunit_module_counter += 1

        # If subunit is empty after Module check, then delete it. 
        if len(newSubunit.modules()) == 0:
            newSubunit.delete()
            #print("***Empty subunit detected: Subunit deleted. ***")
            
        # If subunit is antisense (-), reverse modules in subunit
        if subunitdata[5] == "(-)":
            reverseModuleOrder(newSubunit)
            
    #################### STAGE 3: REORDER SUBUNITS ######################
    # (1) If all subunits are antisense, then reverse the order of every subunit in cluster (e.g. phenalamide1394 or erythromycin54)
    # (2a) Move any subunits with loading modules (AT, ACP only OR CAL only) to the front of the cluster
    # (2b) Move any chain terminating subunits to the end of the cluster
    
    # (1)
    antisenseSubunits = 0
    for finalsubunit in cluster.subunits():
        if finalsubunit.sense == "(-)":
            antisenseSubunits += 1
    if antisenseSubunits == len(cluster.subunits()):
        #print("Reversing subunit order")
        reverseSubunitOrder(cluster)
    
    # Get a list of subunit names
    cluster_subunit_names = []
    for subunit in cluster.subunits():
        cluster_subunit_names.append(subunit.name)
    
    #print(cluster_subunit_names)
    for subunit in cluster.subunits():

        for module in subunit.modules():
            
            # Get a list of module's domains, then check for and move loading modules/subunits
            module_domains = []
            for domain in module.domains():
                module_domains.append(domain.__repr__())
            #print(str(module) + str(module_domains))
            # (2a)
            if module_domains == ['AT','ACP'] or module_domains == ['CAL']:
                #print("SWAP!"+str(module)+str(cluster_subunit_names))
                cluster_subunit_names.remove(subunit.name)
                cluster_subunit_names.insert(0, subunit.name)
                #print("POST SWAP"+str(cluster_subunit_names))
            # (2b) check for and move chain terminating domains to end
            for domain in module.domains():
                if domain.__repr__() == "TE" or domain.__repr__() == "R" or \
                domain.__repr__() == "AOX":
                    cluster_subunit_names.remove(subunit.name)
                    cluster_subunit_names.append(subunit.name)
    
    #print(cluster_subunit_names)
    pks.models.Cluster.reorderSubunits(cluster, cluster_subunit_names)
    cluster.save()
    
    print("[FINAL CLUSTER ARCHITECTURE]")
    for subunit in cluster.subunits():
        print(str(subunit))
        for module in subunit.modules():
            print(str(module) + str(module.domains()))


# If ever needed, Type I subclasses can be accessed at json_data['cluster']['polyketide'].get('synthases')[0]['subclass']
def getPKSNRPSAccessions(filepath):
    filelist = glob.glob(os.path.join(filepath, '*.json'))
    accessions = []
    for filename in filelist:
        with open(filename) as json_file:
            json_data = json.load(json_file)
        try:
            cluster_type = json_data['cluster']['biosyn_class']
            if 'Polyketide' in cluster_type or 'NRP' in cluster_type:
                accession = os.path.basename(filename).strip('.json')
                accessions.append(accession)
        except KeyError:
            pass
        except Exception as e:
            print(e)
            pass
    return accessions