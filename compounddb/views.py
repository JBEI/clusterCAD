from django.http import Http404, HttpResponse
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS, AllChem
from django.utils.http import urlunquote
from django.views.decorators.cache import cache_page
import re
import xml.etree.ElementTree as ET
from io import BytesIO

@cache_page(60 * 120)
def renderer(request, smiles, smiles2=None):
    try:
        # parse smiles input
        smiles = urlunquote(smiles)
        smiles = re.match(r'^(\S{1,10000})', str(smiles)).group(1)
        mol = Chem.MolFromSmiles(smiles)
        if smiles2:
            smiles2 = urlunquote(smiles2)
            smiles2 = re.match(r'^(\S{1,10000})', str(smiles2)).group(1)
    except:
        raise Http404

    if smiles2:
        if ('smarts' in request.GET) and (request.GET['smarts'] == 'True'):
            try:
                template = Chem.MolFromSmarts(smiles2)
            except:
                raise Http404
        else:
            try:
                mol2 = Chem.MolFromSmiles(smiles2)
            except:
                raise Http404
            mcs = rdFMCS.FindMCS([mol, mol2])
            template = Chem.MolFromSmarts(mcs.smartsString)
        highlightAtomLists = mol.GetSubstructMatch(template)
        resultFraction = 100.0 * float(template.GetNumAtoms()) / float(mol.GetNumAtoms())
        legend = "{:.2f}% of atoms".format(resultFraction) 
    else:
        template = Chem.MolFromSmiles('C(=O)[S]')
        highlightAtomLists = None
        legend = ''

    # lock molecule to template
    AllChem.Compute2DCoords(template)
    height = 30 + mol.GetNumAtoms() * 7 
    width = 213 

    align = True 
    if 'align' in request.GET:
        if request.GET['align'] == 'False':
            align = False

    if align:
        try:
            AllChem.GenerateDepictionMatching2DStructure(mol, template, acceptFailure=False)
        except ValueError:
            template = Chem.MolFromSmiles('C(=O)O')
            AllChem.Compute2DCoords(template)
            try:
                AllChem.GenerateDepictionMatching2DStructure(mol, template, acceptFailure=False)
            except ValueError:
                height = 213
    else:
        height = 213

    # draw SVG
    svg = Draw.MolsToGridImage([mol], highlightAtomLists=[highlightAtomLists],
                 legends = [legend], molsPerRow=1,
                 subImgSize=(width, height), useSVG=True)

    # fix XML namespace in RDKit export
    xmlTree = ET.fromstring(svg) 
    xmlTree.attrib['xmlns'] = 'http://www.w3.org/2000/svg'
    xmlTree.attrib['xmlns:svg'] = 'http://www.w3.org/2000/svg'

    # convert XML to string
    fixedSVG = BytesIO()
    ET.ElementTree(xmlTree).write(fixedSVG, encoding='utf-8', xml_declaration=True)

    return HttpResponse(fixedSVG.getvalue(), content_type='image/svg+xml')

