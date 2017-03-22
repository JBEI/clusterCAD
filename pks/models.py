from django.db import models
from model_utils.managers import InheritanceManager
from rdkit import Chem as chem
from rdkit.Chem import AllChem

class Cluster(models.Model):
    # a PKS gene cluster

    genbankAccession = models.CharField(max_length=2000, unique=True) # NCBI accession number
    mibigAccession = models.CharField(max_length=2000, unique=True) # descriptive name
    description = models.TextField() # description of cluster including name of known natural products
    sequence = models.TextField() # nucleotide sequence of cluster 

    def subunits(self):
        # returns a list of all subunits

        return Subunit.objects.filter(cluster=self).order_by('order')

    def architecture(self, products=False):
        # returns the entire structure of the cluster
        # down to individual domains, as a list of lists

        myarchitecture = [[x, x.architecture()] for x in self.subunits()]

        if products:
            try:
                chain = self.computeProduct()
            except:
                chain = False
            i = 0
            for subunit in myarchitecture:
                for module in subunit[1]:
                    if chain:
                        module = module.append(chain[i])
                        i += 1
                    else:
                        module = module.append(False)

        return myarchitecture

    def reorderSubunits(self, newOrder):
        # takes a list of subunit names in order, and reorders
        # them within the gene cluster accordingly

        # confirm that new order is the right length
        # and contains all unique elements
        assert len(newOrder) == len(self.subunits()), 'newOrder has wrong number of subunits for this cluster'
        assert len(newOrder) == len(set(newOrder)), 'newOrder names are not unique'

        subunits = [Subunit.objects.get(name__exact=x, cluster__exact=self) for x in newOrder] 
        for subunit in subunits:
            oldOrder = subunit.order
            modules = subunit.modules()
            subunit.order = None
            subunit.save()
            for module in modules:
                module.subunit = subunit
                module.save()
            subunit.order = oldOrder
            subunit.delete()
        return self.subunits()

    def computeProduct(self):
        chain = []        
        for subunit in self.subunits():
            for module in subunit.modules():
                if len(chain) == 0:
                    chain.append(module.computeProduct())
                else:
                    chain.append(module.computeProduct(chain[-1]))    
        return chain

    def __str__(self):
        return "%s pks gene cluster" % self.description

class Subunit(models.Model):
    # a PKS subunit

    cluster = models.ForeignKey(Cluster)
    order = models.AutoField(primary_key=True)
    genbankAccession = models.CharField(max_length=2000, unique=False) # NCBI accession number
    name = models.CharField(max_length=2000) # name of subunit
    start = models.PositiveIntegerField()
    stop = models.PositiveIntegerField()
    sequence = models.TextField()

    def modules(self):
        return Module.objects.filter(subunit=self).order_by('order')

    def architecture(self):
        return [[x, x.domains()] for x in self.modules()]

    def __str__(self):
        return "%s pks subunit" % self.name

class Module(models.Model):
    # A PKS module within a subunit

    subunit = models.ForeignKey(Subunit)
    order = models.AutoField(primary_key=True)
    loading = models.BooleanField() # Whether or not module is a loading module
    terminal = models.BooleanField() # Whether or not module is a terminal module

    def domains(self):
        return Domain.objects.filter(module=self).select_subclasses().order_by('start')

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

    def computeProduct(self, chain=False):
        domains = {type(domain): domain for domain in self.domains()}
        reactionOrder = [AT, KR, cMT, oMT, DH, ER, TE]
        
        for reaction in reactionOrder:
            if reaction in domains.keys():
                # print("operating with " + str(reaction))
                # if chain:
                #    print(chem.MolToSmiles(chain))
                chain = domains[reaction].operation(chain)

        return chain

    def __str__(self):
        return "pks module %s" % self.order

class Domain(models.Model):
    # this is the parent class for all domains
    # Do not use directly

    module = models.ForeignKey(Module)
    start = models.PositiveIntegerField()
    stop = models.PositiveIntegerField()

    # using InheritanceManager allows us to directly
    # query all Domain subclasses
    objects = InheritanceManager()

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
        ('emal', 'emal'),
        ('mmal', 'mmal'),
        ('mxmal', 'mxmal'),
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
        # if self.module.loading:
        #    if chain:
        #        return chain
        #    else:
        #        return starters[self.substrate]
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

        # if (self.active == True) and (self.module.loading == False):
        if self.active == True:
            assert len(chain.GetSubstructMatches(chem.MolFromSmiles('C(=O)CC(=O)S'),
                       useChirality=True)) == 1, chem.MolToSmiles(chain)
            prod = rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
        else:
            prod = chain

        return prod

    def __str__(self):
        return "type %s" % self.type

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

class KS(Domain):

    def __str__(self):
        return "domain"

class ACP(Domain):

    def __str__(self):
        return "domain"

class PCP(Domain):

    def __str__(self):
        return "domain"

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
