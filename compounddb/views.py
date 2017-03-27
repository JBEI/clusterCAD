from django.http import Http404, HttpResponse
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS, AllChem
from django.utils.http import urlunquote
from django.views.decorators.cache import cache_page
import re
import xml.etree.ElementTree as ET
from io import BytesIO

# @cache_page(60 * 120)
def renderer(request, smiles, smiles2=None):
    try:
        # parse smiles input
        smiles = urlunquote(smiles)
        smiles = re.match(r'^(\S{1,10000})', str(smiles)).group(1)
        mol = Chem.MolFromSmiles(smiles)
        if smiles2:
            smiles2 = urlunquote(smiles2)
            smiles2 = re.match(r'^(\S{1,10000})', str(smiles2)).group(1)
            mol2 = Chem.MolFromSmiles(smiles2)
    except:
        raise Http404

    if smiles2:
        mcs = rdFMCS.FindMCS([mol, mol2])
        template = Chem.MolFromSmarts(mcs.smartsString)
        mols = [mol, mol2]
        highlightAtomLists = [x.GetSubstructMatch(template) for x in mols] 
        molsPerRow = 2
        labels = ['', '']
        width = 213 
    else:
        template = Chem.MolFromSmiles('C(=O)[S]')
        mols = [mol]
        highlightAtomLists = None
        molsPerRow = 1
        labels = ['']
        width = 213 

    # lock molecule to template
    AllChem.Compute2DCoords(template)
    height = 30 + mol.GetNumAtoms() * 7 

    align = True 
    if 'align' in request.GET:
        if request.GET['align'] == 'False':
            align = False

    if align:
        try:
            [AllChem.GenerateDepictionMatching2DStructure(x, template, acceptFailure=False) for x in mols]
        except ValueError:
            template = Chem.MolFromSmiles('C(=O)O')
            AllChem.Compute2DCoords(template)
            try:
                [AllChem.GenerateDepictionMatching2DStructure(x, template, acceptFailure=False) for x in mols]
            except ValueError:
                height = 213
    else:
        height = 213

    # draw SVG
    svg = Draw.MolsToGridImage(mols, highlightAtomLists=highlightAtomLists,
                 legends = labels, molsPerRow=molsPerRow,
                 subImgSize=(width, height), useSVG=True)

    # fix XML namespace in RDKit export
    xmlTree = ET.fromstring(svg) 
    xmlTree.attrib['xmlns'] = 'http://www.w3.org/2000/svg'
    xmlTree.attrib['xmlns:svg'] = 'http://www.w3.org/2000/svg'

    # convert XML to string
    fixedSVG = BytesIO()
    ET.ElementTree(xmlTree).write(fixedSVG, encoding='utf-8', xml_declaration=True)

    return HttpResponse(fixedSVG.getvalue(), content_type='image/svg+xml')

