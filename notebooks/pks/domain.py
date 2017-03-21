from rdkit import Chem as chem
from rdkit.Chem import AllChem, Draw



class Domain(object):
    ''' Abstract base class used to build PKS catalytic domains.
    '''
    def __init__(self):
        pass

# These are the AT substrates supported by antiSMASH:
# {'Malonyl-CoA':'mal',
#  'Methylmalonyl-CoA':'mmal',
#  'Methoxymalonyl-CoA':'mxmal',
#  'Ethylmalonyl-CoA':'emal',
#  'Isobutyryl-CoA':'isobut',
#  '2-Methylbutyryl-CoA':'2metbut',
#  'trans-1,2-CPDA':'trans-1,2-CPDA',
#  'Acetyl-CoA':'Acetyl-CoA',
#  'Benzoyl-_CoA':'benz',
#  'Propionyl-CoA':'prop',
#  '3-Methylbutyryl-CoA':'3metbut',
#  'Ethylmalonyl-CoA':'Ethyl_mal',
#  'CE-Malonyl-CoA':'cemal',
#  '2-Rhyd-Malonyl-CoA':'2Rhydmal',
#  'CHC-CoA':'CHC-CoA',
#  'inactive':'inactive'}

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

class ATL(Domain):
    '''Loading domain initializes polyketide biosynthesis.
    '''    
    def __init__(self, starter_name):
        super(ATL, self).__init__()

        self.starter_name = starter_name
        self.starter_struct = starters[starter_name]

    def operation(self, chain=None):
        '''Optionally takes as input a different starter unit, can be used
           to see what polyketide would be produced if a non-native
           starter unit is used. Otherwise, simply uses native starter.
        '''
        if chain:
            return chain
        else:
            if self.starter_struct:
                return self.starter_struct
            else:
                raise Exception('Please provide a starter unit. The starting unit of this \
                                 loading module is unknown.')

extenders = {'mal': chem.MolFromSmiles('O=C(O)CC(=O)[S]'),
             'mmal': chem.MolFromSmiles('C[C@@H](C(=O)O)C(=O)[S]'),
             'mxmal': chem.MolFromSmiles('CO[C@@H](C(=O)O)C(=O)[S]'),
             'emal': chem.MolFromSmiles('CC[C@@H](C(=O)O)C(=O)[S]')
             }

class KSAT(Domain):
    '''Ketosynthase appends extender unit to nascent acyl chain.
    '''    
    def __init__(self, at_type):
        super(KSAT, self).__init__()

        # Reaction doesn't match stereochemistry at time of query.
        # Needed to change [*:0][C:1] to [C:1], not sure why
        self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[S].'
                                               '[O:4][C:5](=[O:6])[C@@:7][C:8](=[O:9])[S:10]>>'
                                               '[C:1][C:2](=[O:3])[C@:7][C:8](=[O:9])[S:10]'
                                               '.[C:5](=[O:4])(=[O:6])'))
        
        # Assign AT type and cognate extender
        # Currently not allowing arbitrary extender units
        assert at_type in ('mal', 'mmal', 'mxmal', 'emal'), \
                "Invalid AT type '%s' not in (mal, mmal, mxmal, emal)." %(at_type)
        self.at_type = at_type
        # Extender unit specifity is fixed one AT type is assigned
        self.extender = extenders[self.at_type]

    def operation(self, chain):
        assert len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(=O)S'), 
                   useChirality=True)) == 1, chem.MolToSmiles(chain)
    
        prod = self.rxn.RunReactants((chain, self.extender))[0][0]
        chem.SanitizeMol(prod)

        return prod



class KR(Domain):
    ''' Ketoreductase reduces ketone in nascent acyl chain to hydroxyl group.
    '''

    def __init__(self, kr_type, active=True):
        super(KR, self).__init__()

        # Make sure KR type is valid
        assert kr_type in ('A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'U'), \
            "Invalid KR type %s not in ('A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'U')." %(kr_type) 
        # Assign KR type and cognate reaction
        self.kr_type = kr_type
        self.active = active
        # Reaction is fixed once KR type is assigned
        if self.kr_type == 'A1':
            self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C@:2]([O:3])[C@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.kr_type == 'A2':
            self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C@:2]([O:3])[C@@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.kr_type == 'B1':
            self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C@@:2]([O:3])[C@:4]'
                                                   '[C:5](=[O:6])[S:7]'))       
        elif self.kr_type == 'B2':
            self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C@@:2]([O:3])[C@@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.kr_type == 'C1':
            self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.kr_type == 'C2':
            self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C:2](=[O:3])[C@@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        elif self.kr_type == 'C2':
            self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C:2](=[O:3])[C@@:4]'
                                                   '[C:5](=[O:6])[S:7]'))
        else:
            # By first specifying some stereochemistry in the reactants
            # and then explicitly "losing" the stereochemistry in the products
            # we can forget the stereochemistry in our molecule
            self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C@:4]'
                                                   '[C:5](=[O:6])[S:7]>>'
                                                   '[C:1][C:2]([O:3])[C:4]'
                                                   '[C:5](=[O:6])[S:7]'))

    def operation(self, chain):
        assert len(chain.GetSubstructMatches(chem.MolFromSmiles('C(=O)CC(=O)S'), 
                   useChirality=True)) == 1, chem.MolToSmiles(chain)
        if self.active == True:
            prod = self.rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
        else:
            prod = chain
        
        return prod

        
    
class DH(Domain):
    ''' Dehydratase reduces hydroxyl group in nascent acyl chain to methenyl group.
    '''
    def __init__(self, active=True):
        super(DH, self).__init__()

        self.active = active
        self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2]([O:3])[C:4][C:6](=[O:7])[S:8]>>'
                                               '[C:1][CH1:2]=[CH0:4][C:6](=[O:7])[S:8].[O:3]'))

    def operation(self, chain):
#        assert len(chain.GetSubstructMatches(chem.MolFromSmiles('C(O)CC(=O)S'), 
#                   useChirality=True)) == 1, chem.MolToSmiles(chain)

        # Changed assertion to if/else statement so that DH and ER do not 
        # do anything if the KR is inactive
        if len(chain.GetSubstructMatches(chem.MolFromSmiles('C(O)CC(=O)S'), useChirality=True)) == 1 and self.active == True:
            prod = self.rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain


class ER(Domain):
    ''' Enoylreductase reduces methenyl group to methylene group.
    '''
    def __init__(self, active=True):
        super(ER, self).__init__()

        self.active = active
        self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2]=[C:3][C:4](=[O:5])[S:6]>>'
                                               '[C:1][C:2][C@@H1:3][C:4](=[O:5])[S:6]'))

    def operation(self, chain):
#        assert len(chain.GetSubstructMatches(chem.MolFromSmiles('C=CC(=O)S'),
#                   useChirality=True)) == 1, chem.MolToSmiles(chain)
    
        if len(chain.GetSubstructMatches(chem.MolFromSmiles('C=CC(=O)S'), useChirality=True)) == 1 and self.active == True:
            prod = self.rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain



class TE(Domain):
    ''' Thioesterase hyrolyzes and cyclizes nascent acyl chain.
    '''
    def __init__(self, cyclize=False):
        super(TE, self).__init__()

        self.cyclize = cyclize

    def operation(self, chain):
        assert len(chain.GetSubstructMatches(chem.MolFromSmiles('C(=O)S'),
                   useChirality=True)) == 1, chem.MolToSmiles(chain)

        if self.cyclize:
            self.rxn = AllChem.ReactionFromSmarts('([C:1](=[O:2])[S:3].[O:4][C:5][C:6])>>'
                                                  '[C:1](=[O:2])[O:4][C:5][C:6].[S:3]')
        else:
            self.rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])[S:3]>>[C:1](=[O:2])[O].[S:3]')

        # Using index -1 will yield the largest ring
        prod = self.rxn.RunReactants((chain,))[-1][0]
        chem.SanitizeMol(prod)

        return prod



class cMT(Domain):
    def __init__(self, active=True):
        super(cMT, self).__init__()

        self.active = active
        # Note that as implemented, cMT ignores/forgets stereochemistry in alpha position
        self.rxn = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[S:4]>>'
                                               '[C:1](C)[C:2](=[O:3])[S:4]'))

    def operation(self, chain):
        assert len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(=O)S'),
                   )) == 1, chem.MolToSmiles(chain)
    
        if self.active == True:
            prod = self.rxn.RunReactants((chain,))[0][0]
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain



class oMT(Domain):
    def __init__(self, active=True):
        super(oMT, self).__init__()

        self.active = active
        self.rxn1 = AllChem.ReactionFromSmarts(('[C:1][C:2](=[O:3])[C:4][C:5](=[O:6])[S:7]>>'
                                               '[C:1][C:2]([O:3]C)[C:4][C:5](=[O:6])[S:7]'))
        self.rxn2 = AllChem.ReactionFromSmarts(('[C:1][C:2]([O:3])[C:4][C:5](=[O:6])[S:7]>>'
                                               '[C:1][C:2]([O:3]C)[C:4][C:5](=[O:6])[S:7]'))

    def operation(self, chain):
        if self.active == True:
            if len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(=O)CC(=O)S'))) == 1:
                prod = self.rxn1.RunReactants((chain,))[0][0]
            elif len(chain.GetSubstructMatches(chem.MolFromSmiles('CC(O)CC(=O)S'))) == 1: 
                prod = self.rxn2.RunReactants((chain,))[0][0]
            else:
                prod = chain
            chem.SanitizeMol(prod)
            return prod
        else:
            return chain
