from django.shortcuts import render
from django.utils.http import urlunquote, urlquote
from django.contrib import messages
from . import sequencetools
from django.http import Http404
from pks.models import AT, KR

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

    context = {
        'alignments': alignments,
        'aainput': input,
        'evalue': str(evalue),
        'maxHits': str(maxHits),
        'atsubstrates': AT.SUBSTRATE_CHOICES,
        'krtypes': KR.TYPE_CHOICES,
    }
    
    return render(request, 'sequenceresult.html', context)    
    
