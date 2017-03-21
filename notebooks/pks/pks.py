import domain
import numpy as np
from collections import OrderedDict



class Module(object):
    '''Class to represent a PKS module

    # Properties
        domains: dict of dicts. Key is domain name, value is 
                 [dict{'start': start, 'stop': stop}, [specificity notes]],
                 noting that specificity notes are only present for AT and KR.
        loading: bool. Whether or not module is a loading module.
        terminal: bool. Whether or not module is a terminal module.
        operations: list of instances of subclasses of class<Domain>.

    # Methods
        build_operations: Build operations from other attributions
        compute_product: Takes as input a RDKit Mol, and outputs RDKit Mol
                         representing the product of the module
    '''
    def __init__(self, domains, loading=False, terminal=False):
        # Check on validity of domains, should be OrderedDict where keys are modules
        # and values are a list containing a dict of the domain boundaries, which 
        # correspond to the sequence attribute of the parent subunit object
        # and optionally a dict with domain specificity information (only for AT and KR)
        assert isinstance(domains, dict), 'domains must be a dict'
        for d in domains.keys():
            # Currently, cMT and CAL do not do anything
            assert d in ['KS', 'AT', 'KR', 'DH', 'ER', 'ACP',
                         'cMT', 'oMT', 'CAL', 'PCP', 'Thioesterase', 
                         'Heterocyclization', 'AMP-binding', 
                         'Condensation_DCL', 'Condensation_LCL',
                         'PKS_Docking_Nterm', 'PKS_Docking_Cterm'], \
                    'Invalid domain: ' + str(d)
            domain_data = domains[d]
        assert len(set(['ACP', 'PCP']).intersection(set(domains.keys()))) > 0, \
                'Module must contain ACP or PCP:' + str(domains.keys())

        self.domains = domains
        self.loading = loading
        self.terminal = terminal

        self.operations = []
        self.build_operations()

    def build_operations(self):
        "Build operations from attributes."
        # Put together list of chemical operations corresponding to catalytic domains
        operations = []
        if self.loading:
            substratel = self.domains['AT'][1]['Substrate specificity predictions'].split()[0]
            if substratel == 'N/A':
                substratel = 'mal'
            operations.append(domain.ATL(substratel))
            # Stop adding operations if this is a loading module
        else:
            # Inclusion of chemical operations is implemented such that if redundant domains
            # are annotated, then redundant domains will be ignored
            if 'AT' in self.domains.keys():
                 # For now, these checks are done after all Module objects have been instantiated
#                 assert 'KS' in domains.keys(), 'Non-loading module must contain KS: ' + str(domains.keys())
                substrate = self.domains['AT'][1]['Substrate specificity predictions'].split()[0]
                # If the PKS signature prediction is 'N/A', use Minowa prediction
                if substrate == 'N/A':
                    substrate = self.domains['AT'][1]['Substrate specificity predictions'].split()[3]
                operations.append(domain.KSAT(substrate))
            if 'KR' in self.domains.keys():
                stereochemistry =  self.domains['KR'][1]['Predicted KR stereochemistry'].replace('?','U')
                activity = self.domains['KR'][1]['Predicted KR activity']
                if activity == 'active':
                    active = True
                else:
                    active = False
                operations.append(domain.KR(stereochemistry, active))
            if 'cMT' in self.domains.keys():
                operations.append(domain.cMT())
            if 'oMT' in self.domains.keys():
                operations.append(domain.oMT())
            if 'DH' in self.domains.keys():
                operations.append(domain.DH())
            if 'ER' in self.domains.keys():
                operations.append(domain.ER())
            if 'Thioesterase' in self.domains.keys():
                assert self.terminal
                assert self.terminal in set(['Cyclic', 'Linear'])
                if self.terminal == 'Cyclic':
                    operations.append(domain.TE(cyclize=True))
                else:
                    operations.append(domain.TE(cyclize=False))
        self.operations = operations
    
    def compute_product(self, chain):
        '''Function to apply chemical operations of a module to a growing acyl chain.'''
        for operation in self.operations:
            chain = operation.operation(chain)
        return chain



class Subunit(object):
    ''' Class to represent a PKS subunit

    # Properties
        id: str. protein id.  
        name: str. name of subunit.
        description: str. description of the subunit.
        start: int. start of subunit in cluster nucleotide sequence.
        stop: int. end of subunit in cluster nucleotide sequence.
        sequence: str. amino acid sequence corresponding to subunit.
        modules: list of instances of class<Module>.

    # Methods
        compute_product: Takes as input a RDKit Mol, and outputs RDKit Mol
                         representing the product of the subunit.
    '''
    def __init__(self, id, name, description, start, stop, sequence, 
                 modules=[], sspro=None, sspro8=None, accpro=None, accpro20=None):
        self.id = id
        self.name = name
        self.description = description
        self.start = start
        self.stop = stop
        self.sequence = sequence
        self.sspro = sspro
        self.sspro8 = sspro8
        self.accpro = accpro
        self.accpro20 = accpro20
        self.nmodules = len(modules)
        self.modules = modules

    def get_domain_indices(self, module_number, domain_name):
        '''Function that returns indices of a catalytic domain
           contained by a module within amino acid sequence of the subunit.
        '''
        assert module_number < len(self.modules) - 1, \
            'Invalid request, subunit only contains %d modules.' %(len(self.modules))
        module = self.modules[module_number]
        assert domain_name in module.domains.keys(), \
            'Module %d in subunit does not contain catalytic domain %s, \
             options are: %s' %(module_number, domain_name, module.domains.keys())
        return (module[domain_name][0]['start'], module[domain_name][0]['stop'])

    def compute_product(self, chain):
        '''Function to apply chemical operations of a subunit to a growing acyl chain.'''
        for module in self.modules:
            chain = module.compute_product(chain)
        return chain


class Standalone(object):
    ''' Class to represent a PKS standalone enzyme

    # Properties
        id: str. protein id. 
        name: str. name of standalone enzyme.
        description: description of the standalone enzyme.
        start: int. start of standalone enzyme in cluster nucleotide sequence.
        stop: int. end of standalone enzyme in cluster nucleotide sequence.
        sequence: str. amino acid sequence corresponding to standalone enzyme.
    '''
    def __init__(self, id, name, description, start, stop, sequence):
        self.id = id
        self.name = name
        self.description = description
        self.start = start
        self.stop = stop
        self.sequence = sequence



class Cluster(object):
    '''Cluster is defined by ID, name, and description. The subunits that
       comprise the Cluster are provided as list of subunits.

    # Properties
        id: str. NCBI accession number.
        name: str. NCBI accession number.
        description: str. description of cluster.
        sequence: str. nucleotide sequence of cluster.
    '''
    def __init__(self, id, name, description, sequence, 
                 subunits=[], standalones=[]):
        # Information about the cluster
        self.id = id
        self.name = name
        self.description = description
        self.sequence = sequence

        # Input to initilization function expected to be lists of class<Subunit> and class<Standalone>
        # instances respectively; these are converted into OrderedDicts when stored as a cluster object.
        self.subunits = OrderedDict([(subunit.name, subunit) \
                for subunit in subunits]) 
        self.standalones = OrderedDict([(standalone.name, standalone) \
                for standalone in standalones]) 

    def compute_product(self, chain):
        '''Function to apply chemical operations of a cluster to a growing acyl chain.'''
        for subunit in self.subunits.values():
            chain = subunit.compute_product(chain)
        return chain

    def get_starter(self):
        '''Function to return starter unit of cluster.'''
        return list(self.subunits.values())[0].modules[0].operations[0].starter_name

    def update_starter(self, name, structure):
        '''Function to update starter unit of cluster. name is expected to be a string
           representing the starter unit, while structure is expected to be an object
           generated using the chem.MolFromSmiles() function.
        '''
        list(self.subunits.values())[0].modules[0].operations[0].starter_name = name
        list(self.subunits.values())[0].modules[0].operations[0].starter_struct = structure

    def print_subunits(self):
        '''Function to print subunits contained in the cluster.'''
        print(str([subunit.name for subunit in self.subunits.values()]))

    def update_subunit_order(self, order):
        '''Function to update order of subunits in cluster.'''
        # Make sure that all subunits are in input order
        for subunit in self.subunits.keys():
            assert subunit in order, "Missing subunit %s." %(subunit)
        # Recall that order is expected to be a list of the subunit names specifying
        # the new subunit order
        subunit_names = []
        subunit_objects = []
        for subunit_name in order:
            subunit_names.append(subunit_name)
            subunit_objects.append(self.subunits[subunit_name])
        updated_subunits = OrderedDict([pair for pair in zip(subunit_names, subunit_objects)])           # Update old starter module
        updated_subunits[list(self.subunits.keys())[0]].modules[0].loading = False
        updated_subunits[list(self.subunits.keys())[0]].modules[0].build_operations()
        # Update new starter module
        updated_subunits[order[0]].modules[0].loading = True
        updated_subunits[order[0]].modules[0].build_operations()
        self.subunits = updated_subunits

    def pop_subunit(self, name):
        del self.subunits[name]

    def toggle_cyclization(self):
        last_module = self.subunits[list(self.subunits.keys())[-1]].modules[-1]
        if 'Thioesterase' in last_module.domains.keys():
            old_cyclization = last_module.operations[-1].cyclize
            print('Old cyclization specification was %s' %(old_cyclization))
            if old_cyclization == False:
                new_cyclization = True
            else: 
                new_cyclization = False
            print('New cyclization specification is %s' %(new_cyclization))
            self.subunits[list(self.subunits.keys())[-1]].modules[-1].operations[-1].cyclize = new_cyclization
        else:
            print('Last subunit does not have a TE!')
            print(last_module.domains.keys())
            print(last_module.operations)

    def print_domains(self):
        '''Function to print catalytic domains contained in the cluster.'''
        for subunit in list(self.subunits.values()):
            print(subunit.name)
            for module in subunit.modules:
                print(str(list(module.domains.keys())))


#    def get_module(self, module):
#        '''Function to get a module contained in the cluster. Domain should be queried by
#           providing the number corresponding to the module of interest within the cluster. 
#        '''
#        # This gets the subunits as a list
#        subunits = list(self.subunits.values())
#        # Counts of the modules in each subunit, make sure requested module is contained in cluster
#        module_counts = np.cumsum(subunit.nmodules for subunit in self.subunits)
#        assert module < module_counts[-1], 'Cluster only contains %d modules.' %(module_counts[-1])
#        # Figure out which subunit contains the module of interest
#        subunit_index = 0
#        while module < module_counts[subunit_index]
#            subunit_index += 1
#        # Get module
#        module = subunits[subunit_index].modules[module - module_counts[subunit_index-1]]
#        return module

