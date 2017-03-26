from django.shortcuts import render
from django.utils.http import urlunquote
from django.contrib import messages
from compounddb.models import Compound
from pks.models import Module
from rdkit import Chem as chem

def search(request):
    if request.method != 'POST':
        if 'smiles' in request.GET:
            smiles = urlunquote(request.GET['smiles'])
        else:
            smiles = ''
        context = {'smiles': smiles}
        return render(request, 'search.html', context)
    elif 'smiles' in request.POST:
        try:
            mol = smilesToMol(request.POST['smiles'])
        except:
            messages.error(request, 'Error: Invalid SMILES string!')
            return render(request, 'search.html')
    elif 'sdf' in request.POST:
        try:
            mol = chem.MolFromMolBlock(request.POST['sdf'])
            # assert mol != None
            if not mol:
                raise ValueError
        except ValueError:
            messages.error(request, 'Error: Invalid structure!')
            messages.error(request, request.POST['sdf'])
            return render(request, 'search.html')

    compoundHits = Compound.atomPairSearch(chem.MolToSmiles(mol), minSim=0.0)

    moduleHits = [] 
    for similarity, compound in compoundHits:
        modules = list(Module.objects.filter(product=compound))
        if len(modules) == 0:
            continue 
        hitDict = {
            'similarity': similarity,
            'smiles': compound.smiles,
            'modules': modules
        }
        moduleHits.append(hitDict)

    context = {'moduleHits': moduleHits}
    
    return render(request, 'result.html', context)    
    

def smilesToMol(smiles):
    mol = chem.MolFromSmiles(smiles)
    chem.SanitizeMol(mol)
    return mol 
