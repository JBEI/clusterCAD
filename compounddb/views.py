from django.http import Http404, HttpResponse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from django.utils.http import urlunquote
from django.views.decorators.cache import cache_page
import re
import xml.etree.ElementTree as ET
from io import BytesIO

# @cache_page(60 * 120)
def renderer(request, smiles, acpDisplace=0):
    try:
        # parse smiles input
        smiles = urlunquote(smiles)
        smiles = re.match(r'^(\S{1,10000})', str(smiles)).group(1)
        mol = Chem.MolFromSmiles(smiles)
        acpDisplace = float(acpDisplace)

    except:
        raise Http404

    # lock molecule to template
    template = Chem.MolFromSmiles('C(=O)[S]')
    AllChem.Compute2DCoords(template)
    height = 30 + mol.GetNumAtoms() * 7 

    onACP = False
    
    align = True 
    if 'align' in request.GET:
        if request.GET['align'] == 'False':
            align = False

    if align:
        try:
            AllChem.GenerateDepictionMatching2DStructure(mol, template, acceptFailure=False)
            onACP = True
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
    width = 213 
    svg = Draw.MolsToGridImage([mol],
                 legends = [''], molsPerRow=1,
                 subImgSize=(width, height), useSVG=True)

    # fix XML namespace in RDKit export
    xmlTree = ET.fromstring(svg) 
    xmlTree.attrib['xmlns'] = 'http://www.w3.org/2000/svg'
    xmlTree.attrib['xmlns:svg'] = 'http://www.w3.org/2000/svg'

    # center sulfur atom 
    # by making S predictable, it can be lined up with the ACP
    onACP = False # disable this feature for now
    rect = [child for child in xmlTree[0] if child.tag == 'rect'][0]
    rect.attrib['style'] = 'opacity:0.0;stroke:none'
    if onACP:
        s = [child for child in xmlTree[0] if child.tag == 'text' and child[0].text == 'S'][0]
        sX = float(s.attrib['x'])
        sY = float(s.attrib['y'])
        newX = acpDisplace - sX
        newY = 10 - sY
        if newX < 0:
            newX = 0
        xmlTree[0].attrib['transform'] = 'translate(%s,%s)' % (newX, newY,)
        rect.attrib['width'] = str(width + newX) 
        xmlTree.attrib['width'] = '%spx' % rect.attrib['width']

    # convert XML to string
    fixedSVG = BytesIO()
    ET.ElementTree(xmlTree).write(fixedSVG, encoding='utf-8', xml_declaration=True)

    return HttpResponse(fixedSVG.getvalue(), content_type='image/svg+xml')

