from django.http import Http404, HttpResponse
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS, AllChem
from django.utils.http import urlunquote
from django.views.decorators.cache import cache_page
import re
import xml.etree.ElementTree as ET
from io import BytesIO

# @cache_page(60 * 120)
def renderer(request, smiles, width=243):
    try:
        # parse smiles input
        smiles = urlunquote(smiles)
        smiles = re.match(r'^(\S{1,10000})', str(smiles)).group(1)
        mol = Chem.MolFromSmiles(smiles)
        width=int(width)
        assert 0 < width < 2000
    except:
        raise Http404

    template = Chem.MolFromSmiles('C(=O)[S]')
    highlightAtomLists = None
    legend = ''

    # lock molecule to template
    AllChem.Compute2DCoords(template)
    height = int(width * ((30.0 + mol.GetNumAtoms() * 7.0) / 243.0))

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
                height = width 
    else:
        height = width

    # draw SVG
    svg = Draw.MolsToGridImage([mol], highlightAtomLists=[highlightAtomLists],
                 molsPerRow=1,
                 subImgSize=(width, height), useSVG=True)

    # fix XML namespace in RDKit export
    xmlTree = ET.fromstring(svg) 
    xmlTree.attrib['xmlns'] = 'http://www.w3.org/2000/svg'
    xmlTree.attrib['xmlns:svg'] = 'http://www.w3.org/2000/svg'

    # convert XML to string
    fixedSVG = BytesIO()
    ET.ElementTree(xmlTree).write(fixedSVG, encoding='utf-8', xml_declaration=True)

    return HttpResponse(fixedSVG.getvalue(), content_type='image/svg+xml')

# @cache_page(60 * 120)
def mcsrenderer(request, smiles1, smiles2, mcsSmarts=False, width=852, align=True, chiral=False):
    # parse input
    try:
        smiles1 = urlunquote(smiles1)
        smiles1 = re.match(r'^(\S{1,10000})', str(smiles1)).group(1)
        mol1 = Chem.MolFromSmiles(smiles1)
        smiles2 = urlunquote(smiles2)
        smiles2 = re.match(r'^(\S{1,10000})', str(smiles2)).group(1)
        mol2 = Chem.MolFromSmiles(smiles2)
        width=int(width)
        assert 0 < width < 2000
        align=int(align)
        assert 0 <= align <= 1
        chiral=int(chiral)
        assert 0 <= chiral <= 1
    except:
        raise Http404

    # find MCS if not provided
    if mcsSmarts and mcsSmarts != '':
        mcsSmarts = urlunquote(mcsSmarts)
        mcsSmarts = re.match(r'^(\S{1,10000})', str(mcsSmarts)).group(1)
        template = Chem.MolFromSmarts(mcsSmarts)
    else:
        if chiral:
            mcs = rdFMCS.FindMCS([mol1, mol2], matchChiralTag=True)
        else:
            mcs = rdFMCS.FindMCS([mol1, mol2])
        template = Chem.MolFromSmarts(mcs.smartsString)
    highlightAtomList1 = mol1.GetSubstructMatch(template)
    highlightAtomList2 = mol2.GetSubstructMatch(template)
    
    if align:
        # generate coords for first molecule 
        AllChem.Compute2DCoords(mol1)
        atoms = [atom.GetIdx() for atom in mol1.GetAtoms()]
        nonMCSatoms = list(filter(lambda x: x not in highlightAtomList1, atoms))
        nonMCSatoms.sort(reverse=True)
        rwmol = Chem.RWMol(mol1)
        for atom in nonMCSatoms:
            rwmol.RemoveAtom(atom)
        templateMol = rwmol.GetMol()

        # apply coords to second molecule
        AllChem.GenerateDepictionMatching2DStructure(mol2, templateMol, acceptFailure=True)

    # draw SVG
    width = int(float(width) / 2)
    height = int(float(width) * (350.0 / 426.0))
    svg = Draw.MolsToGridImage([mol1, mol2], highlightAtomLists=[highlightAtomList1, highlightAtomList2],
                 molsPerRow=2,
                 subImgSize=(width, height), useSVG=True)

    # fix XML namespace in RDKit export
    xmlTree = ET.fromstring(svg) 
    xmlTree.attrib['xmlns'] = 'http://www.w3.org/2000/svg'
    xmlTree.attrib['xmlns:svg'] = 'http://www.w3.org/2000/svg'

    # convert XML to string
    fixedSVG = BytesIO()
    ET.ElementTree(xmlTree).write(fixedSVG, encoding='utf-8', xml_declaration=True)

    return HttpResponse(fixedSVG.getvalue(), content_type='image/svg+xml')
