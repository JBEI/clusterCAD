from django.shortcuts import render
from django.utils.http import urlunquote, urlquote
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
            mol = chem.MolFromMolBlock(request.POST['sdf'], strictParsing=False)
            # assert mol != None
            if not mol:
                raise ValueError
        except ValueError:
            messages.error(request, 'Error: Invalid structure!')
            return render(request, 'search.html')

    querySmiles = chem.MolToSmiles(mol, isomericSmiles=True)

    assert 'cutoff' in request.POST
    minSim = float(request.POST['cutoff'])
    assert 0.0 <= minSim <= 1.0

    assert 'maxCompounds' in request.POST
    maxHits = int(request.POST['maxCompounds'])
    assert 0 < maxHits <= 1000 

    compoundHits = Compound.atomPairSearch(querySmiles, minSim=minSim, maxHits=maxHits)

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

    if len(moduleHits) == 0:
        messages.error(request, 'No hits - please refine query!')
        return render(request, 'search.html')

    context = {
        'querySmiles': urlquote(querySmiles),
        'moduleHits': moduleHits,
        'minSim': minSim,
        'maxHits': maxHits
    }
    
    return render(request, 'result.html', context)    
    

def smilesToMol(smiles):
    mol = chem.MolFromSmiles(smiles)
    chem.SanitizeMol(mol)
    return mol 
