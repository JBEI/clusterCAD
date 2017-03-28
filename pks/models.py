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

    # Methods
        subunits: Returns subunits in cluster.
        architecture: Returns structure of cluster as list of lists.
        reorderSubunits: Reorders subunits in cluster.
        computeProduct: Compute list of intermediates for each module in cluster. 
        setActive: Sets module KR, DH, ER, oMT, or cMT domain active boolean.
        setSubstrate: Sets AT substrate specificity of module.
        setStereochemistry: Sets KR stereochemistry of module.
        setCyclization: Sets cyclization of cluster TE.
        clusterDict: Generates OrderedDict representation of cluster.
        clusterJSON: Generates JSON file representation of cluster and save to file path.
        correctCluster: Corrects errors in cluster using JSON file. 
                        Template JSON file should be generated using clusterJSON.
    '''
    mibigAccession = models.CharField(max_length=100, primary_key=True)
    genbankAccession = models.CharField(max_length=100)
    description = models.TextField()
    sequence = models.TextField()
    knownProductSmiles = models.TextField()
    knownProductMCS = models.TextField()
    knownProductSource = models.TextField(default='unknown')

    def subunits(self):
        return Subunit.objects.filter(cluster=self).order_by('order')

    def architecture(self):
        # Returns the entire structure of the cluster
        return [[x, x.architecture()] for x in self.subunits()]

    def reorderSubunits(self, newOrder):
        '''Takes as input a list of subunit names and reordereds subunits within
           the gene cluster.
        '''
        for subunit in self.subunits():
            assert subunit.name in newOrder, 'Missing subunit %s.' %(subunit)
        subunits = [Subunit.objects.get(name__exact=x, cluster__exact=self) for x in newOrder] 
        assert len(newOrder) == len(subunits), 'Non-existant subunits provided in new order.'

        # Reset loading bool of first module in oldOrder
        oldLoading = self.subunits()[0].modules()[0]
        oldLoading.loading = False
        oldLoading.setLoading()
        oldLoading.save()

        # Update subunit ordering
        subunitCounter = 1
        moduleCounter = 1
        for subunit in subunits:
            subunit.order = subunitCounter
            print(subunit.order)
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
        newLoading = self.subunits()[0].modules()[0]
        newLoading.loading = True
        newLoading.setLoading()
        newLoading.save()

        # compute new final products after reordering
        self.computeProduct()

        return self.subunits()

    def computeProduct(self, computeMCS=True):
        chain = []        
        for subunit in self.subunits():
            for module in subunit.modules():
                if len(chain) == 0:
                    chain.append(module.computeProduct())
                else:
                    chain.append(module.computeProduct(chain[-1]))    

        if computeMCS:
                mcs = rdFMCS.FindMCS([chem.MolFromSmiles(self.knownProductSmiles), chain[-1]])
                self.knownProductMCS = mcs.smartsString
                self.save()
        return chain

    def setActive(self, sub, mod, dom, active):
        assert dom in ['KR', 'DH', 'ER']
        assert isinstance(active, bool)
        domain = Domain.objects.filter(module__subunit__cluster=self, 
                                       module__subunit__name=sub).select_subclasses( \
                                       getattr(sys.modules[__name__], dom)).order_by( \
                                       'start')[mod]
        domain.active = active
        domain.save()

    def setSubstrate(self, sub, mod, update):
        assert update in ['mal', 'mmal', 'mxmal', 'emal', 'cemal',
                          'prop', 'isobut', '2metbut', 'trans-1,2-CPDA', 'CHC-CoA']
        d = AT.objects.filter(module__subunit__cluster=self, 
                              module__subunit__name=sub).order_by('start')[mod]
        d.substrate = update
        d.save()

    def setStereochemistry(self, sub, mod, update):
        assert update in ['A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'U']
        d = KR.objects.filter(module__subunit__cluster=self,
                              module__subunit__name=sub).order_by('start')[mod]
        d.type = update
        d.save()

    def setCyclization(self, cyclic):
        d = TE.objects.get(module__subunit__cluster=self)
        d.cyclic = cyclic
        d.save()

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
                    elif dname == 'KR':
                        ddict.update({'active': domain.active,
                                      'type': domain.type})
                    elif dname in ['DH', 'ER', 'cMT', 'oMT']:
                        ddict.update({'active': domain.active})
                    elif dname == 'TE':
                        ddict.update({'cyclic': domain.cyclic})
                    if dname in ['AT', 'KR', 'DH', 'ER', 'cMT', 'oMT', 'TE']:
                        mdict.update({dname: ddict})
                sdict.update({imodule: mdict})
            adict.update({sname: sdict})
        ret.update({'architecture': adict})
        return ret

    def clusterJSON(self, path):
        '''Function that writes a JSON file containing the changeable parameters
            describing a PKS cluster.
        '''
        with open(os.path.join(path, self.mibigAccession + '.json'), 'w') as f:
            f.write(json.dumps(self.clusterDict(), indent=4))

    def correctCluster(self, filepath):
        '''Function that takes as input a JSON file representing the changeable parameters
           describing a PKS cluster and updates the database entry to reflect the information
           contained in the JSON file.
        '''
        corr = json.loads(open(filepath).read(), object_pairs_hook=OrderedDict)
        assert corr['mibigAccession'] == self.mibigAccession
        assert corr['genbankAccession'] == self.genbankAccession
        # Reorder subunits if necessary
        newOrder = [str(x) for x in corr['architecture'].keys()]
        self.reorderSubunits(newOrder)
        # Change domain properties if necessary
        for s,sdict in corr['architecture'].items():
            for m,mdict in sdict.items():
                m = int(m) # Key is expected to be integer
                for d,ddict in mdict.items():
                    if d == 'AT':
                        self.setSubstrate(s, m, ddict['substrate'])
                    elif d == 'KR':
                        self.setStereochemistry(s, m, ddict['type'])
                        self.setActive(s, m, d, ddict['active'])
                    elif d in ['DH', 'ER']:
                        self.setActive(s, m, d, ddict['active'])
                    else:
                        assert d == 'TE'
                        cyclic = ddict['cyclic']
                        self.setCyclization(cyclic)
    
    def __str__(self):
        return "%s gene cluster" % self.description

class Subunit(models.Model):
    ''' Class representing a PKS subunit.

    # Properties
        cluster: class<Cluster>. cluster containing subunit.
        order: int. order of subunit within cluster starting with 1. Auto-set on self.save()
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
    cluster = models.ForeignKey(Cluster)
    order = models.IntegerField()
    genbankAccession = models.CharField(max_length=2000, unique=False)
    name = models.CharField(max_length=2000)
    start = models.PositiveIntegerField()
    stop = models.PositiveIntegerField()
    sequence = models.TextField()

    def modules(self):
        return Module.objects.filter(subunit=self).order_by('order')
    
    def architecture(self):
        return [[x, x.domains()] for x in self.modules()]

    def getNucleotideSequence(self):
        return self.cluster.sequence[self.start:self.stop]

    def getAminoAcidSequence(self):
        return self.sequence

    def __str__(self):
        return "%s" % self.name

@receiver(pre_save, sender=Subunit)
def setSubunitOrder(sender, instance, **kwargs):
    # sets the subunit order
    if not instance.order:
        subunitCount = sender.objects.filter(cluster=instance.cluster).count()
        instance.order = subunitCount + 1 

class Module(models.Model):
    '''Class representing a PKS module.

    # Properties
        subunit: class<Subunit>. subunit containing module.
        order: int. order of module within cluster starting with 1. Auto-set on self.save()
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
    subunit = models.ForeignKey(Subunit)
    order = models.IntegerField()
    loading = models.BooleanField() # Whether or not module is a loading module
    terminal = models.BooleanField() # Whether or not module is a terminal module
    product = models.ForeignKey(Compound, on_delete=models.SET_NULL, default=None, blank=True, null=True) # small molecule product structure

    def domains(self):
        return Domain.objects.filter(module=self).select_subclasses().order_by('start')

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
        if 'KS' in domainDict.keys():
            start = domainDict['KS'][0]['start']
            stop = domainDict['KS'][0]['stop']
            newDomain = KS(module=self, start=start, stop=stop)
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

        if 'cMT' in domainDict.keys():
            start = domainDict['cMT'][0]['start']
            stop = domainDict['cMT'][0]['stop']
            newDomain = cMT(module=self, start=start, stop=stop)
            newDomain.save()

        if 'oMT' in domainDict.keys():
            start = domainDict['oMT'][0]['start']
            stop = domainDict['oMT'][0]['stop']
            newDomain = oMT(module=self, start=start, stop=stop)
            newDomain.save()

        if 'ACP' in domainDict.keys():
            start = domainDict['ACP'][0]['start']
            stop = domainDict['ACP'][0]['stop']
            newDomain = ACP(module=self, start=start, stop=stop)
            newDomain.save()

        if 'PCP' in domainDict.keys():
            start = domainDict['PCP'][0]['start']
            stop = domainDict['PCP'][0]['stop']
            newDomain = PCP(module=self, start=start, stop=stop)
            newDomain.save()

        if 'Thioesterase' in domainDict.keys():
            start = domainDict['Thioesterase'][0]['start']
            stop = domainDict['Thioesterase'][0]['stop']
            newDomain = TE(module=self, start=start, stop=stop, cyclic=cyclic)
            newDomain.save()
        
        self.setLoading()

    def computeProduct(self, chain=False):
        if self.product:
            return self.product.mol()
        domains = {type(domain): domain for domain in self.domains()}
        reactionOrder = [AT, KR, cMT, oMT, DH, ER, TE]
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
        if Module.objects.filter(product=self.product).exclude(id=self.id).count() == 0:
            self.product.delete()
        self.product = None

    def __str__(self):
        return "%s" % str(self.order)

@receiver(pre_delete, sender=Module)
def deleteModuleProduct(sender, instance, **kwargs):
    instance.deleteProduct()

@receiver(pre_save, sender=Module)
def setModuleOrder(sender, instance, **kwargs):
    # sets the module order based on current count
    if not instance.order:
        moduleCount = sender.objects.filter(subunit__cluster=instance.subunit.cluster).count()
        instance.order = moduleCount + 1 

class Domain(models.Model):
    ''' Abstract base class used to build PKS catalytic domains.

    # Properties
        module: class<Module>. module containing domain.
        start: int. start of domain in subunit amino acid sequence.
        stop: int. end of domain in subunit amino acid sequence.

    # Methods
        getNucleotideSequence: Returns nucleotide sequence of domain.
        getAminoAcidSequence: Returns amino acid sequence of domain.
    '''
    module = models.ForeignKey(Module)
    start = models.PositiveIntegerField()
    stop = models.PositiveIntegerField()

    # Using InheritanceManager allows us to directly
    # query all Domain subclasses
    objects = InheritanceManager()

    def getNucleotideSequence(self):
        sequence = self.module.subunit.getNucleotideSequence()
        return sequence[self.start*3:self.stop*3]

    def getAminoAcidSequence(self):
        sequence = self.module.subunit.getAminoAcidSequence()
        return sequence[self.start:self.stop]

# dict of supported starter units
starters = {'mal': chem.MolFromSmiles('CC(=O)[S]'),
            'mmal': chem.MolFromSmiles('CCC(=O)[S]'),
            'mxmal': chem.MolFromSmiles('COCC(=O)[S]'),
            'cemal': chem.MolFromSmiles('CC(=O)[S]'),
            'Acetyl-CoA': chem.MolFromSmiles('CC(=O)[S]'),
            'prop': chem.MolFromSmiles('CCC(=O)[S]'),
            'isobut': chem.MolFromSmiles('CC(C)C(=O)[S]'),
            '2metbut': chem.MolFromSmiles('CCC(C)C(=O)[S]'),
            'CHC-CoA': chem.MolFromSmiles('C1CCCCC1C(=O)[S]'),
            'trans-1,2-CPDA': chem.MolFromSmiles('C1[C@@H](C(=O)O)CC[C@@H]1C(=O)[S]'),
            'N/A': None
           }

extenders = {'mal': chem.MolFromSmiles('O=C(O)CC(=O)[S]'),
             'mmal': chem.MolFromSmiles('C[C@@H](C(=O)O)C(=O)[S]'),
             'mxmal': chem.MolFromSmiles('CO[C@@H](C(=O)O)C(=O)[S]'),
             'emal': chem.MolFromSmiles('CC[C@@H](C(=O)O)C(=O)[S]')
             }

class AT(Domain):
    SUBSTRATE_CHOICES = (
        ('mal', 'mal'),
        ('mmal', 'mmal'),
        ('mxmal', 'mxmal'),
        ('emal', 'emal'),
        ('cemal', 'cemal'),
        ('Acetyl-CoA', 'Acetyl-CoA'),
        ('prop', 'prop'),
        ('isobut', 'isobut'),
        ('2metbut', '2metbut'),
        ('CHC-CoA', 'CHC-CoA'),
        ('trans-1,2-CPDA', 'trans-1,2-CPDA'),
        ('N/A', 'N/A'),
    )
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
            rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[S].'
                                              '[O:4][C:5](=[O:6])[C@@:7][C:8](=[O:9])[S:10]>>'
                                              '[C:1][C:2](=[O:3])[C@:7][C:8](=[O:9])[S:10]'
                                              '.[C:5](=[O:4])(=[O:6])'))
            assert KS in [type(domain) for domain in self.module.domains()]
            assert len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(=O)S'),
                       useChirality=True)) == 1, chem.MolToSmiles(chain)
            prod = rxn.RunReactants((chain, extenders[self.substrate]))[0][0]
            chem.SanitizeMol(prod)
            return prod

    def __str__(self):
        return "substrate %s, loading %s" % (self.substrate, self.module.loading)

    def __repr__(self):
        return("AT")

class KR(Domain):
    active = models.BooleanField()
    TYPE_CHOICES = (
        ('A1', 'A1'),
        ('A2', 'A2'),
        ('B1', 'B1'),
        ('B2', 'B2'),
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
            rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C@:2]([O:3])[C@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'A2':
            rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C@:2]([O:3])[C@@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'B1':
            rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C@@:2]([O:3])[C@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'B2':
            rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C@@:2]([O:3])[C@@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'C1':
            rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.type == 'C2':
            rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C:2](=[O:3])[C@@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        else:
            # By first specifying some stereochemistry in the reactants
            # and then explicitly "losing" the stereochemistry in the products
            # we can forget the stereochemistry in our molecule
            rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C@:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C:2]([O:3])[C:4]'
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
        return "type %s, active %s" % (self.type, self.active)

    def __repr__(self):
        return("KR")

class DH(Domain):
    active = models.BooleanField(default=True)

    def operation(self, chain):
        rxn = AllChem.ReactionFromSmarts(('[C:1][C:2]([O:3])[C:4][C:6](=[O:7])[S:8]>>'
                                          '[C:1][CH1:2]=[CH0:4][C:6](=[O:7])[S:8].[O:3]'))
        # Changed assertion to if/else statement so that DH and ER do not
        # do anything if the KR is inactive
        if len(chain.GetSubstructMatches(chem.MolFromSmiles('C(O)CC(=O)S'), useChirality=True)) == 1 and self.active == True:
            prod = rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain

    def __str__(self):
        return "active %s" % self.active

    def __repr__(self):
        return("DH")

class ER(Domain):
    active = models.BooleanField(default=True)

    def operation(self, chain):
        rxn = AllChem.ReactionFromSmarts(('[C:1][C:2]=[C:3][C:4](=[O:5])[S:6]>>'
                                          '[C:1][C:2][C@@H1:3][C:4](=[O:5])[S:6]'))
        if len(chain.GetSubstructMatches(chem.MolFromSmiles('C=CC(=O)S'), useChirality=True)) == 1 and self.active == True:
            prod = rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain

    def __str__(self):
        return "active %s" % self.active

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
        return "active %s" % self.active

    def __repr__(self):
        return("cMT")

class oMT(Domain):
    active = models.BooleanField(default=True)

    def operation(self, chain):
        if self.active == True:
            if len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(=O)CC(=O)S'))) == 1:
                rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4][C:5](=[O:6])[S:7]>>'
                                                  '[C:1][C:2]([O:3]C)[C:4][C:5](=[O:6])[S:7]'))
                prod = self.rxn.RunReactants((chain,))[0][0]
            elif len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(O)CC(=O)S'))) == 1:
                rxn = AllChem.ReactionFromSmarts(('[C:1][C:2]([O:3])[C:4][C:5](=[O:6])[S:7]>>'
                                                  '[C:1][C:2]([O:3]C)[C:4][C:5](=[O:6])[S:7]'))
                prod = self.rxn.RunReactants((chain,))[0][0]
            else:
                prod = chain
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain

    def __str__(self):
        return "active %s" % self.active

    def __repr__(self):
        return("oMT")

class TE(Domain):
    cyclic = models.BooleanField()

    def operation(self, chain):
        assert len(chain.GetSubstructMatches(chem.MolFromSmiles('C(=O)S'),
                   useChirality=True)) == 1, chem.MolToSmiles(chain)
        if self.cyclic:
            rxn = AllChem.ReactionFromSmarts('([C:1](=[O:2])[S:3].[O:4][C:5][C:6])>>'
                                                  '[C:1](=[O:2])[O:4][C:5][C:6].[S:3]')
        else:
            rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])[S:3]>>[C:1](=[O:2])[O].[S:3]')

        # Using index -1 will yield the largest ring
        prod = rxn.RunReactants((chain,))[-1][0]
        chem.SanitizeMol(prod)

        return prod

    def __str__(self):
        return "cyclic %s" % self.cyclic

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
    cluster = models.ForeignKey(Cluster)
    order = models.PositiveSmallIntegerField()
    name = models.CharField(max_length=2000) # name of enzyme 
    start = models.PositiveIntegerField()
    stop = models.PositiveIntegerField()
    sequence = models.TextField()

# API example code:

'''
import pks.models

# create a gene cluster
myCluster = pks.models.Cluster(genbankAccession = "myGB", genbankVersion="1", mibigAccession="myMB", mibigVersion="1", description="test cluster", sequence="ACTG")
myCluster.save()

# create and save subunits
# note: ordering is automatically determined by the order
# in which they are saved
subunit1 = pks.models.Subunit(cluster=myCluster, name="s1", start=1, stop=3, sequence="AAC")
subunit1.save()
subunit2 = pks.models.Subunit(cluster=myCluster, name="s2", start=4, stop=4, sequence="AAC")
subunit2.save()

# create and save modules
# note: ordering is automatically determined by the order
# in which they are saved
myModule1 = pks.models.Module(subunit=subunit1, loading=True, terminal=False)
myModule1.save()
myModule2 = pks.models.Module(subunit=subunit1, loading=False, terminal=True)
myModule2.save()

# create and save domains
myATL = pks.models.ATL(module=myModule1, start=1, stop=20, starter='mmal')
myATL.save()
myER = pks.models.ER(module=myModule1, start=21, stop=22, active=True)
myER.save()

# get cluster from database
myCluster = pks.models.Cluster.objects.all()[0]
myCluster.architecture() # show full structure of cluster
myCluster.subunits()
mySubunit = myCluster.subunits()[0]
mySubunit.modules()
myModule1 = mySubunit.modules()[0]
myModule1.domains()
myDomain = myModule1.domains()[0]

# query domain properties
myDomain.start
myDomain.stop
myDomain.starter
'''
