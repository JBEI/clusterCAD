from django.shortcuts import render
from django.utils.http import urlunquote, urlquote
from django.contrib import messages
from . import sequencetools
from django.http import Http404
from pks.models import AT, KR, DH, ER, cMT, oMT, TE

def search(request):
    if request.method != 'POST':
        if 'aainput' in request.GET:
            aainput = urlunquote(request.GET['aainput'])
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
            alignments = sequencetools.blast(query=input, evalue=evalue, max_target_seqs=maxHits)
        except:
            messages.error(request, 'Error: Invalid query!')
            return render(request, 'sequencesearch.html')
    else:
        raise Http404

    if len(alignments) == 0:
        messages.error(request, 'No hits - please refine query!')
        return render(request, 'sequencesearch.html')

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
        'aainput': input,
        'evalue': str(evalue),
        'maxHits': str(maxHits),
        'atsubstrates': atsubstrates,
        'krtypes': krtypes,
        'boolDomains': boolDomains,
        'tetypes': tetypes,
    }
    
    return render(request, 'sequenceresult.html', context)    
    
