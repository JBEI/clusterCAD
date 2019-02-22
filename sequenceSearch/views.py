import re
from django.shortcuts import render
from django.utils.http import urlunquote, urlquote
from django.contrib import messages
from . import sequencetools
from django.http import Http404
from model_utils.managers import InheritanceManager
from django.core.cache import cache

from pks.models import AT, KR, DH, ER, cMT, oMT, TE, Subunit, Domain

def search(request):
    if request.method != 'POST':
        if 'aainput' in request.GET:
            aainput = urlunquote(request.GET['aainput'])
            messages.success(request, 'Imported AA sequence from previous page')
        elif 'subunit' in request.GET:
            subunitid = int(request.GET['subunit'])
            assert 0 <= subunitid
            subunit = Subunit.objects.get(id=subunitid)
            aainput = subunit.sequence
            messages.success(request, 'Imported sequence for PKS cluster %s subunit %s' % (subunit.cluster.description, subunit.name))
        else:
            aainput = ''
        context = {'aainput': aainput}
        return render(request, 'sequencesearch.html', context)
    elif 'aainput' in request.POST:
        try:
            input = request.POST['aainput']
            maxHits = int(request.POST['maxHits'])
            assert 1 <= maxHits <= 10000
            evalue = float(request.POST['evalue'])
            assert 0.0 <= evalue <= 10.0
            showAllDomains = int(request.POST['showAllDomains'])
            assert 0 <= showAllDomains <= 1
            input = input.strip()
            if len(re.findall('>.*?\n', input)) >= 2:
                messages.error(request, "Error: Multiple queries detected, please remove until only one query is present")
                return render(request, 'sequencesearch.html')
            input = re.sub('^>.*?\n', '', input)
            input = re.sub('\n', '', input)
            input = re.sub('\s', '', input)
            if len(input) > 50000:
                messages.error(request, 'Error: max query size is 50,000 residues')
                return render(request, 'sequencesearch.html')

            # use alignment cache if it exists
            alignments = cache.get((input, evalue, maxHits))
            if not alignments:
                alignments = sequencetools.blast(query=input, evalue=evalue, max_target_seqs=maxHits)
                if alignments == 0:
                    messages.error(request, 'Error: Query not found')
                    return render(request, 'sequencesearch.html')
                cache.set((input, evalue, maxHits), alignments, 60 * 60 * 24 * 7) # cache for one week
                
        except ValueError:
            messages.error(request, 'Error: Invalid query!')
            return render(request, 'sequencesearch.html')
    else:
        raise Http404

    if len(alignments) == 0:
        messages.error(request, 'No hits - please refine query!')
        return render(request, 'sequencesearch.html')

    # if showAllDomains = 1, we will replace the contents of 
    # each domain list with all domains in it's module
    if showAllDomains:
        for alignment in alignments:
            for hsp in alignment['hsps']:
                for module in hsp['modules']:
                    module['domains'] = list(Domain.objects.filter(module=module['module']).select_subclasses().order_by('start'))

    # get domain options to display in UI
    domains = [domain for alignment in alignments for hsp in alignment['hsps'] for module in hsp['modules'] for domain in module['domains']]
    ats = list(filter(lambda d: isinstance(d, AT), domains))
    atsubstrates = list(set([at.substrate for at in ats]))
    krs = list(filter(lambda d: isinstance(d, KR), domains))
    krtypes = list(set([kr.type for kr in krs]))
    boolDomains = []
    for domain in (DH, ER, cMT, oMT):
        thesedomains = list(filter(lambda d: isinstance(d, domain), domains))
        typelist = list(set([str(d) for d in thesedomains]))
        if len(typelist) > 0:
            boolDomains.append((domain.__name__, typelist))
    tes = list(filter(lambda d: isinstance(d, TE), domains))
    tetypes = list(set([str(te) for te in tes]))

    context = {
        'alignments': alignments,
        'queryResidues': len(input),
        'evalue': str(evalue),
        'maxHits': str(maxHits),
        'showAllDomains': showAllDomains,
        'atsubstrates': atsubstrates,
        'krtypes': krtypes,
        'boolDomains': boolDomains,
        'tetypes': tetypes,
    }
    
    return render(request, 'sequenceresult.html', context)    
    
