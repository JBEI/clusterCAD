import sys
from django.db import models
from model_utils.managers import InheritanceManager
from rdkit import Chem as chem
from rdkit.Chem import AllChem, rdFMCS
import os
import json
from collections import OrderedDict
from compounddb.models import Compound
from django.db.models.signals import pre_save, pre_delete
from django.dispatch import receiver
from  django.core.validators import validate_comma_separated_integer_list
class Cluster(models.Model):
    '''Cluster is defined by GenBank and MIBiG accession numbers.

    # Properties
        mibigAccession: str. MIBiG accession number.
        genbankAccession: str. GenBank accession number.
        description: str. known natural product of cluster.
        sequence: str. nucleotide sequence of cluster.
        knownProductSmiles: str. Known final product from literature or external databases, in smiles format.
        knownProductMCS: str. SMARTS string of MCS comparing knownProductSmiles to the predicted final product. Auto-generated when running self.computeProduct().
        knownProductSource: str. Where the knownProduct structure was obtained from.
        reviewed: bool. True if the final chemical structure has been manually reviewed for correctness. This is usually set to True only if a corrections file is provided.

    # Methods
        subunits: Returns subunits in cluster.
        architecture: Returns structure of cluster as list of lists.
        reorderSubunits: Reorders subunits in cluster.
        computeProduct: Compute list of intermediates for each module in cluster.
        setIterations: Sets module iterations parameter specificying number of iterations.
        setActive: Sets module KR, DH, ER, oMT, or cMT domain active boolean.
        setSubstrate: Sets AT substrate specificity of module.
        setStereochemistry: Sets KR stereochemistry of module.
        setCyclization: Sets cyclization of cluster TE.
        clusterDict: Generates OrderedDict representation of cluster.
        clusterJSON: Generates JSON file representation of cluster and save to file path.
        correctCluster: Corrects errors in cluster using JSON file.
                        Template JSON file should be generated using clusterJSON.
    '''
    #id = models.autofield(primary_key=True)
    mibigAccession = models.CharField(max_length=100)
    genbankAccession = models.CharField(max_length=100)
    description = models.TextField()
    sequence = models.TextField()
    hasKnownProduct = models.BooleanField(default=False)
    knownProductSmiles = models.TextField(default=False)
    knownProductMCS = models.TextField()
    knownProductSource = models.TextField(default='unknown')
    reviewed = models.BooleanField(default=False)

    def subunits(self):
        if self.isOrdered():
            return Subunit.objects.filter(cluster=self).order_by('order')
        return Subunit.objects.filter(cluster=self)

    def architecture(self):
        # Returns the entire structure of the cluster
        return [[x, x.architecture()] for x in self.subunits()]

    def isOrdered(self):
        return all(subunit.isOrdered() for subunit in Subunit.objects.filter(cluster=self))

    def computeProduct(self, computeMCS=True, recompute=False):
        if not self.isOrdered():
            return []
            #raise Exception('Chemical Product prediction not available for pks with unknown subunit ordering.'))
        chain = []
        for subunit in self.subunits():
            for module in subunit.modules():
                if recompute:
                    module.deleteProduct()
                if len(chain) == 0:
                    chain.append(module.computeProduct())
                else:
                    chain.append(module.computeProduct(chain[-1]))

        if computeMCS:
                mcs = rdFMCS.FindMCS([chem.MolFromSmiles(self.knownProductSmiles), chain[-1]])
                self.knownProductMCS = mcs.smartsString
                self.save()
        return chain

    def reorderSubunits(self, newOrder):
        '''Takes as input a list of subunit names and reorders subunits within
           the gene cluster.
        '''
        for subunit in self.subunits():
            assert subunit.name in newOrder, 'Missing subunit %s.' %(subunit)
        subunits = [Subunit.objects.get(cluster__exact=self, name__exact=x) for x in newOrder]
        assert len(newOrder) == len(subunits), 'Non-existant subunits provided in new order.'

        # Reset loading bool of first module in oldOrder
        if len(self.subunits()[0].modules()) > 0:
            oldLoading = self.subunits()[0].modules()[0]
            oldLoading.loading = False
            oldLoading.setLoading()
            oldLoading.save()

        # Update subunit ordering
        subunitCounter = 0
        moduleCounter = 0
        for subunit in subunits:
            subunit.order = subunitCounter
            subunit.save()
            subunitCounter += 1

            # update module ordering
            modules = subunit.modules()
            for module in modules:
                module.order = moduleCounter
                module.deleteProduct()
                module.save()
                moduleCounter += 1

        # Reset loading bool of first module in newOrder
        newLoading = list(filter(lambda subunit: subunit.order == 0, self.subunits()))[0].modules()[0]
        newLoading.loading = True
        newLoading.setLoading()
        newLoading.save()

        return self.subunits()

    def setIterations(self, sub, mod, iterations):
        module = Module.objects.filter(subunit__cluster=self,
                                       subunit__name=sub).order_by('order')[mod]
        module.iterations = iterations
        module.save()

    def setActive(self, sub, mod, dom, active):
        assert dom in ['KR', 'DH', 'ER', 'oMT', 'cMT']
        assert isinstance(active, bool)
        module = Module.objects.filter(subunit__cluster=self,
                                       subunit__name=sub).order_by('order')[mod]
        domains = Domain.objects.filter(module=module).select_subclasses()
        domain = list(filter(lambda x: repr(x) == dom, list(domains)))[0]
        domain.active = active
        domain.save()

    def setSubstrate(self, sub, mod, update):
        assert update in list(set(list(starters.keys()) + list(extenders.keys()))), update
        module = Module.objects.filter(subunit__cluster=self,
                                       subunit__name=sub).order_by('order')[mod]
        try:
            domain = AT.objects.get(module=module)
        except AT.DoesNotExist:
            domain = CAL.objects.get(module=module)
        domain.substrate = update
        domain.save()

    def setStereochemistry(self, sub, mod, update):
        assert update in [x[0] for x in KR.TYPE_CHOICES]
        module = Module.objects.filter(subunit__cluster=self,
                                       subunit__name=sub).order_by('order')[mod]
        domain = KR.objects.get(module=module)
        domain.type = update
        domain.save()

    def setCyclization(self, cyclic, ring=0):
        assert isinstance(cyclic, bool)
        assert isinstance(ring, int)
        domain = TE.objects.get(module__subunit__cluster=self)
        domain.cyclic = cyclic
        if cyclic:
            domain.ring = ring
        domain.save()

    def clusterDict(self):
        '''Function that generates OrderedDict representation of changeable parameters
           describing a PKS cluster.
        '''
        ret = OrderedDict()
        ret.update({'description': self.description})
        ret.update({'genbankAccession': self.genbankAccession})
        ret.update({'mibigAccession': self.mibigAccession})
        adict = OrderedDict()
        for subunit in self.subunits():
            sname = subunit.name
            sdict = OrderedDict()
            for imodule,module in enumerate(subunit.modules()):
                mdict = OrderedDict()
                for domain in module.domains():
                    dname = repr(domain)
                    ddict = OrderedDict()
                    if dname == 'AT':
                        ddict.update({'substrate': domain.substrate})
                    elif dname == 'CAL':
                        ddict.update({'substrate': domain.substrate})
                    elif dname == 'KR':
                        ddict.update({'active': domain.active,
                                      'type': domain.type})
                    elif dname in ['DH', 'ER', 'cMT', 'oMT']:
                        ddict.update({'active': domain.active})
                    elif dname == 'TE':
                        ddict.update({'cyclic': domain.cyclic,
                                      'ring': domain.ring})
                    if dname in ['AT', 'CAL', 'KR', 'DH', 'ER', 'cMT', 'oMT', 'TE']:
                        mdict.update({dname: ddict})
                sdict.update({imodule: {'domains': mdict,
                                        'iterations': module.iterations}})
            adict.update({sname: sdict})
        ret.update({'architecture': adict})
        return ret

    def clusterJSON(self, path):
        '''Function that writes a JSON file containing the changeable parameters
            describing a PKS cluster.
        '''
        with open(os.path.join(path, self.id + '.json'), 'w') as f:
            f.write(json.dumps(self.clusterDict(), indent=4))

    def correctCluster(self, filepath):
        '''Function that takes as input a JSON file representing the changeable parameters
           describing a PKS cluster and updates the database entry to reflect the information
           contained in the JSON file.
        '''
        corr = json.loads(open(filepath).read(), object_pairs_hook=OrderedDict)
        assert corr['mibigAccession'] == self.mibigAccession
        assert corr['genbankAccession'] == self.genbankAccession
        # Delete subunits if necessary
        for s in Subunit.objects.filter(cluster=self):
            if s.name not in corr['architecture'].keys():
                Subunit.objects.get(cluster=self, name=s).delete()
        # Change domain properties if necessary
        for s,sdict in corr['architecture'].items():
            for m,mdict_iters in sdict.items():
                m = int(m) # Key is expected to be integer
                mdict = mdict_iters['domains']
                iters = mdict_iters['iterations']
                self.setIterations(s, m, iters)
                for d,ddict in mdict.items():
                    if d == 'AT':
                        self.setSubstrate(s, m, ddict['substrate'])
                    elif d == 'CAL':
                        self.setSubstrate(s, m, ddict['substrate'])
                    elif d == 'KR':
                        self.setStereochemistry(s, m, ddict['type'])
                        self.setActive(s, m, d, ddict['active'])
                    elif d in ['DH', 'ER', 'cMT', 'oMT']:
                        self.setActive(s, m, d, ddict['active'])
                    else:
                        assert d == 'TE'
                        cyclic = ddict['cyclic']
                        ring = int(ddict['ring'])
                        self.setCyclization(cyclic, ring)
        # Reorder subunits if necessary
        # Must be done last since this is when products are recomputed
        newOrder = [str(x) for x in corr['architecture'].keys()]
        self.reorderSubunits(newOrder)
        self.computeProduct(recompute=True)

    def __str__(self):
        if len(self.description) == 0:
            return 'unnamed gene cluster'#?
        return "%s gene cluster" % self.description

class Subunit(models.Model):
    ''' Class representing a PKS subunit.

    # Properties
        cluster: class<Cluster>. cluster containing subunit.
        order: int. order of subunit within cluster starting with 0. Auto-set on self.save(). Set to -1 if not known.
        genbankAccession: str. GenBank accession number.
        name: str. name of subunit.
        start: int. start of subunit in cluster nucleotide sequence.
        stop: int. end of subunit in cluster nucleotide sequence.
        sequence: str. amino acid sequence corresponding to subunit.

    # Methods
        modules: Returns modules in subunit.
        domains: Returns domains in subunit.
        architecture: Returns architecture of subunit.
        getNucleotideSequence: Returns nucleotide sequence of subunit.
        getAminoAcidSequence: Returns amino acid sequence of subunit.
    '''
    cluster = models.ForeignKey(Cluster, on_delete=models.CASCADE)
    order = models.IntegerField()
    genbankAccession = models.CharField(max_length=2000, unique=False)
    name = models.CharField(max_length=2000)
    start = models.PositiveIntegerField()
    stop = models.PositiveIntegerField()
    sequence = models.TextField()
    acc = models.TextField()
    acc20 = models.CharField(validators = [validate_comma_separated_integer_list], max_length=1000000)
    accPlotFile = models.ImageField(upload_to='accplots')
    ss = models.TextField()
    ss8 = models.TextField()
    ssPlotFile = models.ImageField(upload_to='ssplots')

    def modules(self):
        return Module.objects.filter(subunit=self).select_subclasses().order_by('order')

    def architecture(self):
        return [[x, x.domains()] for x in self.modules()]

    def getNucleotideSequence(self):
        return self.cluster.sequence[(self.start - 1):self.stop]

    def getAminoAcidSequence(self):
        return self.sequence

    def isOrdered(self):
        return self.order != -1

    def __str__(self):
        return "%s" % self.name

@receiver(pre_save, sender=Subunit)
def setSubunitOrder(sender, instance, **kwargs):
    # sets the subunit order
    if not isinstance(instance.order, int):
        subunitCount = sender.objects.filter(cluster=instance.cluster).count()
        instance.order = subunitCount


class DomainContainer(models.Model):
    """
    Abstract base class for anything that contains a domain. Offers the domains()
    method, which returns a list of domains.
    """
    def domains(self):
        return Domain.objects.filter(container=self).select_subclasses().order_by('start')
    class Meta:
        abstract = True


class Standalone(models.Model):
    # a standalone PKS enzyme within a gene cluster
    '''
    A standalone PKS enzyme with specific domains.
    # Properties
        cluster: class<Cluster>. cluster containing subunit.
        genbankAccession: str. GenBank accession number.
        name: str. name of standaloen.
        start: int. start of subunit in cluster nucleotide sequence.
        stop: int. end of subunit in cluster nucleotide sequence.
        sequence: str. amino acid sequence corresponding to subunit.
    # Methods
        architecture: Returns architecture of subunit.
        getNucleotideSequence: Returns nucleotide sequence of subunit.
        getAminoAcidSequence: Returns amino acid sequence of subunit.
    '''
    cluster = models.ForeignKey(Cluster, on_delete=models.CASCADE)
    name = models.CharField(max_length=2000) # name of enzyme
    start = models.PositiveIntegerField()
    stop = models.PositiveIntegerField()
    sequence = models.TextField()
    genbankAccession = models.CharField(max_length=100)

class DomainContainingStandalone(Standalone, DomainContainer):
    """
    A standalone PKS enzyme with specific domains.
    """
    pass

class Module(DomainContainer):
    '''Class representing a PKS module.

    # Properties
        subunit: class<Subunit>. subunit containing module.
        order: int. order of module within cluster starting with 0. Auto-set on self.save()
        loading: bool. Whether or not module is a loading module.
        terminal: bool. Whether or not module is a terminal module.
        product: class<Compound>. Product small molecule structure.

    # Methods
        domains: Returns domains in subunit.
        setLoading: Resets activity of reductive casette based on whether module is loading.
        buildDomains: Build class<Domain> objects using dict as input
        computeProduct: Compute product of module given chain.
        deleteProduct: Reset product to Null, and properly delete it from database.
    '''
    subunit = models.ForeignKey(Subunit, on_delete=models.CASCADE)
    order = models.IntegerField()
    loading = models.BooleanField() # Whether or not module is a loading module
    terminal = models.BooleanField() # Whether or not module is a terminal module
    product = models.ForeignKey(Compound, on_delete=models.SET_NULL, default=None, blank=True, null=True) # small molecule product structure
    iterations = models.PositiveIntegerField(default=1)
    objects = InheritanceManager()

    def setLoading(self):
        domains = Domain.objects.filter(module=self).select_subclasses(KR, DH, ER, cMT, oMT)
        if self.loading:
            for d in domains:
                d.active = False
                d.save()
        else:
            for d in domains:
                d.active=True
                d.save()

    def buildDomains(self, domainDict, cyclic=False):
        if 'CAL' in domainDict.keys():
            start = domainDict['CAL'][0]['start']
            stop = domainDict['CAL'][0]['stop']
            substrate = domainDict['CAL'][1]['Substrate specificity predictions'].split()[0]
            if substrate == 'N/A':
                substrate = domainDict['CAL'][1]['Substrate specificity predictions'].split()[3]
            if substrate == 'N/A':
                substrate = 'mal'
            newDomain = CAL(module=self, start=start, stop=stop, substrate=substrate)
            newDomain.save()

        if 'AT' in domainDict.keys():
            start = domainDict['AT'][0]['start']
            stop = domainDict['AT'][0]['stop']
            substrate = domainDict['AT'][1]['Substrate specificity predictions'].split()[0]
            if substrate == 'N/A':
                substrate = domainDict['AT'][1]['Substrate specificity predictions'].split()[3]
            if substrate == 'N/A':
                substrate = 'mal'
            newDomain = AT(module=self, start=start, stop=stop, substrate=substrate)
            newDomain.save()

        if 'KR' in domainDict.keys():
            start = domainDict['KR'][0]['start']
            stop = domainDict['KR'][0]['stop']
            type =  domainDict['KR'][1]['Predicted KR stereochemistry']
            if type == '?':
                type = 'U'
            activity = domainDict['KR'][1]['Predicted KR activity']
            if activity == 'active':
                active = True
            else:
                active = False
            newDomain = KR(module=self, start=start, stop=stop, active=active, type=type)
            newDomain.save()

        for domainType in ['KS', 'DH', 'ER', 'cMT', 'oMT', 'ACP', 'PCP']:
            if domainType in domainDict.keys():
                start = domainDict[domainType][0]['start']
                stop = domainDict[domainType][0]['stop']
                if domainType in ['DH', 'ER', 'cMT', 'oMT']:
                    newDomain = getattr(sys.modules[__name__], domainType)(module=self, start=start, stop=stop, active=True)
                else:
                    newDomain = getattr(sys.modules[__name__], domainType)(module=self, start=start, stop=stop)
                newDomain.save()

        if 'Thioesterase' in domainDict.keys():
            start = domainDict['Thioesterase'][0]['start']
            stop = domainDict['Thioesterase'][0]['stop']
            newDomain = TE(module=self, start=start, stop=stop, cyclic=cyclic, ring=0)
            newDomain.save()

        self.setLoading()

    def computeProduct(self, chain=False):
        if self.product:
            return self.product.mol()
        domains = {type(domain): domain for domain in self.domains()}
        reactionOrder = [CAL, AT, KR, cMT, oMT, DH, ER, TE]
        for iteration in range(self.iterations):
            for reaction in reactionOrder:
                if reaction in domains.keys():
                    chain = domains[reaction].operation(chain)
        # Save this modules product in the database
        thisProduct = Compound(smiles = chem.MolToSmiles(chain, isomericSmiles=True))
        thisProduct.save()
        self.product = thisProduct
        self.save()
        return chain

    def deleteProduct(self):
        # set self.product to none, and delete the compound itself
        # if this is the only module
        if Module.objects.filter(product=self.product).exclude(id=self.id).count() == 0 and self.product != None:
            self.product.delete()
        self.product = None

    def loadingStr(self):
        if self.loading:
            return 'loading'
        else:
            return 'non-loading'

    def __str__(self):
        return "%s" % str(self.order)

@receiver(pre_delete, sender=Module)
def deleteModuleProduct(sender, instance, **kwargs):
    instance.deleteProduct()

@receiver(pre_save, sender=Module)
def setModuleOrder(sender, instance, **kwargs):
    # sets the module order based on current count
    if not isinstance(instance.order, int):
        moduleCount = sender.objects.filter(subunit__cluster=instance.subunit.cluster).count()
        instance.order = moduleCount

class Domain(models.Model):
    ''' Abstract base class used to build PKS catalytic domains.

    # Properties
        container: class<DomainContainer>. module containing domain.
        start: int. start of domain in subunit amino acid sequence.
        stop: int. end of domain in subunit amino acid sequence.

    # Methods
        getNucleotideSequence: Returns nucleotide sequence of domain.
        getAminoAcidSequence: Returns amino acid sequence of domain.
    '''
    container = models.ForeignKey(DomainContainer, on_delete=models.CASCADE, default=None, null=True)
    start = models.PositiveIntegerField()
    stop = models.PositiveIntegerField()

    # Using InheritanceManager allows us to directly
    # query all Domain subclasses
    objects = InheritanceManager()

    def getAminoAcidSequence(self):
        sequence = self.module.subunit.getAminoAcidSequence()
        return sequence[(self.start - 1):self.stop]

    class Meta:
        abstract = True

def activityString(domain):
    if domain.active:
        return 'active'
    else:
        return 'inactive'

# dict of supported starter units
starters = {'mal': chem.MolFromSmiles('CC(=O)[S]'),
            'mmal': chem.MolFromSmiles('CCC(=O)[S]'),
            'butmal': chem.MolFromSmiles('CCCC(=O)[S]'),
            'mxmal': chem.MolFromSmiles('COCC(=O)[S]'),
            'mxmal_ACP': chem.MolFromSmiles('COCC(=O)[S]'),
            'cemal': chem.MolFromSmiles('CC(=O)[S]'),
            'Acetyl-CoA': chem.MolFromSmiles('CC(=O)[S]'),
            'prop': chem.MolFromSmiles('CCC(=O)[S]'),
            'isobut': chem.MolFromSmiles('CC(C)C(=O)[S]'),
            '2metbut': chem.MolFromSmiles('CCC(C)C(=O)[S]'),
            'CHC-CoA': chem.MolFromSmiles('C1CCCCC1C(=O)[S]'),
            'trans-1,2-CPDA': chem.MolFromSmiles('C1CC[C@@H](C(=O)O)[C@@H]1C(=O)[S]'),
            'cyclopentene': chem.MolFromSmiles('C1(=O)C(=CCC1)C(=O)[S]'),
            'pyr': chem.MolFromSmiles('[S]C(C1=CC=CN1)=O'),
            'cin': chem.MolFromSmiles('O=C([S])/C=C/C1=CC=CC=C1'),
            'AHBA': chem.MolFromSmiles('[S]C(C1=CC(O)=CC(N)=C1)=O'),
            'isoval': chem.MolFromSmiles('CC(CC([S])=O)C'),
            'PABA': chem.MolFromSmiles('NC1=CC=C(C([S])=O)C=C1'),
            'DCP': chem.MolFromSmiles('[S]C(C1=CNC(Cl)=C1Cl)=O'),
            'guan': chem.MolFromSmiles('NC(NCC([S])=O)=[NH2+]'),
            'mthz': chem.MolFromSmiles('CC1=NC(C([S])=O)=CS1'),
            'DHCH': chem.MolFromSmiles('O[C@H]1[C@H](O)CCC(C([S])=O)C1'),
            'DHCHene': chem.MolFromSmiles('O[C@H]1[C@H](O)CC=C(C([S])=O)C1'),
            'plac': chem.MolFromSmiles('O=C([S])CC1=CC=CC=C1'),
            'benz': chem.MolFromSmiles('[S]C(C1=CC=CC=C1)=O'),
            'PNBA': chem.MolFromSmiles('[S]C(C1=CC=C([N+]([O-])=O)C=C1)=O'),
            'ema': chem.MolFromSmiles('[S]C([C@@H](CC)C(N)=O)=O'),
            '3measp': chem.MolFromSmiles('[S]C([C@@H](C)CNC([C@@H](N)C)=O)=O'),
            # the following structures still need confirmation
            'AAnon': chem.MolFromSmiles('[S]C(CC(NC([C@H](C)N)=O)CCCCCC)=O'),
            'a3abut': chem.MolFromSmiles('[S]C(C[C@@H](NC(C)=O)C)=O'),
            'shikimic_acid': chem.MolFromSmiles('C1[C@H]([C@@H]([C@@H](C=C1C(=O)[S])O)O)O'),
            'fatty_acid': chem.MolFromSmiles('CCCC(=O)[S]'),
            'NH2': chem.MolFromSmiles('NCC(=O)[S]'),
            'N/A': None
           }

extenders = {'mal': chem.MolFromSmiles('O=C(O)CC(=O)[S]'),
             'mmal': chem.MolFromSmiles('C[C@@H](C(=O)O)C(=O)[S]'),
             'mxmal': chem.MolFromSmiles('CO[C@@H](C(=O)O)C(=O)[S]'),
             'mxmal_ACP': chem.MolFromSmiles('CO[C@@H](C(=O)O)C(=O)[S]'),
             'emal': chem.MolFromSmiles('CC[C@@H](C(=O)O)C(=O)[S]'),
             'allylmal': chem.MolFromSmiles('C=CC[C@@H](C(=O)O)C(=O)[S]'),
             'butmal': chem.MolFromSmiles('CCCC[C@@H](C(=O)O)C(=O)[S]'),
             'hmal': chem.MolFromSmiles('OC(C(C([S])=O)O)=O'),
             'isobutmal': chem.MolFromSmiles('CC(C)[C@@H](C([S])=O)C(O)=O'),
             'D-isobutmal': chem.MolFromSmiles('CC(C)[C@H](C([S])=O)C(O)=O'),
             'DCP': chem.MolFromSmiles('ClC1=C(Cl)NC=C1CCCCC(C(O)=O)C([S])=O'),
             'hexmal': chem.MolFromSmiles('CCCCCC[C@@H](C(=O)O)C(=O)[S]')
             }

class CAL(Domain):
    SUBSTRATE_CHOICES = tuple(((k, k) for k in list(set(list(extenders.keys()) + list(starters.keys())))))

    substrate = models.CharField(
        max_length=20,
        choices=SUBSTRATE_CHOICES,
        default='mal',
        blank=False,
    )

    def operation(self, chain):
        assert self.module.loading == True
        return starters[self.substrate]

    def __str__(self):
        if self.module.iterations > 1:
            return "substrate %s, %s, iterations %d" % (self.substrate, self.module.loadingStr(), self.module.iterations)
        else:
            return "substrate %s, %s" % (self.substrate, self.module.loadingStr())

    def __repr__(self):
        return("CAL")

class AT(Domain):
    SUBSTRATE_CHOICES = tuple(((k, k) for k in list(set(list(extenders.keys()) + list(starters.keys())))))

    substrate = models.CharField(
        max_length=20,
        choices=SUBSTRATE_CHOICES,
        default='mal',
        blank=False,
    )

    def operation(self, chain):
        if not chain:
            return starters[self.substrate]
        else:
            rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[S].'
                                              '[O:4][C:5](=[O:6])[C@@:7][C:8](=[O:9])[S:10]>>'
                                              '[C:2](=[O:3])[C@:7][C:8](=[O:9])[S:10]'
                                              '.[C:5](=[O:4])(=[O:6])'))
            assert KS in [type(domain) for domain in self.module.domains()]
            assert len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(=O)S'),
                       useChirality=True)) == 1, chem.MolToSmiles(chain)
            prod = rxn.RunReactants((chain, extenders[self.substrate]))[0][0]
            chem.SanitizeMol(prod)
            return prod

    def __str__(self):
        if self.module.iterations > 1:
            return "substrate %s, %s, iterations %d" % (self.substrate, self.module.loadingStr(), self.module.iterations)
        else:
            return "substrate %s, %s" % (self.substrate, self.module.loadingStr())

    def __repr__(self):
        return("AT")

class KR(Domain):
    active = models.BooleanField()
    TYPE_CHOICES = (
        ('A1', 'A1'),
        ('A2', 'A2'),
        ('A', 'A'),
        ('B1', 'B1'),
        ('B2', 'B2'),
        ('B', 'B'),
        ('C1', 'C1'),
        ('C2', 'C2'),
        ('U', 'U'),
    )
    type = models.CharField(
        max_length=2,
        choices=TYPE_CHOICES,
        default=None,
        blank=False,
    )

    def operation(self, chain):
        # Reaction is fixed once KR type is assigned
        if self.type == 'A1':
            rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C@:2]([O:3])[C@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'A2':
            rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C@:2]([O:3])[C@@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'A':
            rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C@:2]([O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'B1':
            rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C@@:2]([O:3])[C@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'B2':
            rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C@@:2]([O:3])[C@@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'B':
            rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C@@:2]([O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'C1':
            rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'C2':
            rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:2](=[O:3])[C@@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        else:
            # By first specifying some stereochemistry in the reactants
            # and then explicitly "losing" the stereochemistry in the products
            # we can forget the stereochemistry in our molecule
            rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[C@:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:2]([O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        if self.active == True:
            assert len(chain.GetSubstructMatches(chem.MolFromSmiles('C(=O)CC(=O)S'),
                       useChirality=True)) == 1, chem.MolToSmiles(chain)
            prod = rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
        else:
            prod = chain
        return prod

    def __str__(self):
        return "type %s, %s" % (self.type, activityString(self))

    def __repr__(self):
        return("KR")

class DH(Domain):
    active = models.BooleanField(default=True)

    def operation(self, chain):
        rxn = AllChem.ReactionFromSmarts(('[C:2]([O:3])[C:4][C:6](=[O:7])[S:8]>>'
                                          '[CH1:2]=[CH0:4][C:6](=[O:7])[S:8].[O:3]'))
        # Changed assertion to if/else statement so that DH and ER do not
        # do anything if the KR is inactive
        if len(chain.GetSubstructMatches(chem.MolFromSmiles('C(O)CC(=O)S'), useChirality=True)) == 1 and self.active == True:
            prod = rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain

    def __str__(self):
        return activityString(self)

    def __repr__(self):
        return("DH")

class ER(Domain):
    active = models.BooleanField(default=True)

    def operation(self, chain):
        rxn = AllChem.ReactionFromSmarts(('[C:2]=[C:3][C:4](=[O:5])[S:6]>>'
                                          '[C:2][C@@H1:3][C:4](=[O:5])[S:6]'))
        if len(chain.GetSubstructMatches(chem.MolFromSmiles('C=CC(=O)S'), useChirality=True)) == 1 and self.active == True:
            prod = rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain

    def __str__(self):
        return activityString(self)

    def __repr__(self):
        return("ER")

class cMT(Domain):
    active = models.BooleanField(default=True)

    def operation(self, chain):
        rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[S:4]>>'
                                          '[C:1](C)[C:2](=[O:3])[S:4]'))
        assert len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(=O)S'),
                   )) == 1, chem.MolToSmiles(chain)
        if self.active == True:
            prod = rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain

    def __str__(self):
        return activityString(self)

    def __repr__(self):
        return("cMT")

class oMT(Domain):
    active = models.BooleanField(default=True)

    def operation(self, chain):
        if self.active == True:
            if len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(=O)CC(=O)S'))) == 1:
                rxn = AllChem.ReactionFromSmarts(('[C:2](=[O:3])[C:4][C:5](=[O:6])[S:7]>>'
                                                  '[C:2]([O:3]C)[C:4][C:5](=[O:6])[S:7]'))
                prod = rxn.RunReactants((chain,))[0][0]
            elif len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(O)CC(=O)S'))) == 1:
                rxn = AllChem.ReactionFromSmarts(('[C:2]([O:3])[C:4][C:5](=[O:6])[S:7]>>'
                                                  '[C:2]([O:3]C)[C:4][C:5](=[O:6])[S:7]'))
                prod = rxn.RunReactants((chain,))[0][0]
            else:
                prod = chain
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain

    def __str__(self):
        return activityString(self)

    def __repr__(self):
        return("oMT")

class TE(Domain):
    cyclic = models.BooleanField()
    ring = models.IntegerField(default=0)

    def operation(self, chain):
        assert len(chain.GetSubstructMatches(chem.MolFromSmiles('C(=O)S'),
                   useChirality=True)) == 1, chem.MolToSmiles(chain)

        index = -1
        if self.cyclic:
            rxn = AllChem.ReactionFromSmarts('([C:1](=[O:2])[S:3].[O:4][C:5][C:6])>>'
                                                  '[C:1](=[O:2])[O:4][C:5][C:6].[S:3]')
            index -= self.ring
        else:
            rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])[S:3]>>[C:1](=[O:2])[O].[S:3]')

        # Using index -1 will yield the largest ring
        prod = rxn.RunReactants((chain,))[index][0]
        chem.SanitizeMol(prod)

        return prod

    def __str__(self):
        if self.cyclic:
            return 'cyclic'
        else:
            return 'non-cyclic'

    def __repr__(self):
        return("TE")

class KS(Domain):

    def __str__(self):
        return "domain"

    def __repr__(self):
        return("KS")

class ACP(Domain):

    def __str__(self):
        return "domain"

    def __repr__(self):
        return("ACP")

class PCP(Domain):

    def __str__(self):
        return "domain"

    def __repr__(self):
        return("PCP")

class Standalone(models.Model):
    # a standalone PKS enzyme within a gene cluster
    '''
    A standalone PKS enzyme with specific domains.
    # Properties
        cluster: class<Cluster>. cluster containing subunit.
        genbankAccession: str. GenBank accession number.
        name: str. name of standaloen.
        start: int. start of subunit in cluster nucleotide sequence.
        stop: int. end of subunit in cluster nucleotide sequence.
        sequence: str. amino acid sequence corresponding to subunit.
    # Methods
        architecture: Returns architecture of subunit.
        getNucleotideSequence: Returns nucleotide sequence of subunit.
        getAminoAcidSequence: Returns amino acid sequence of subunit.
    '''
    cluster = models.ForeignKey(Cluster, on_delete=models.CASCADE)
    name = models.CharField(max_length=2000) # name of enzyme
    start = models.PositiveIntegerField()
    stop = models.PositiveIntegerField()
    sequence = models.TextField()
    genbankAccession = models.CharField(max_length=100)

