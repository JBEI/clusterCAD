import re
from django.shortcuts import render
from django.utils.http import urlunquote, urlquote
from django.contrib import messages
from sequenceSearch.tasks import blast
from django.http import Http404, JsonResponse
from model_utils.managers import InheritanceManager
from pks.models import *
import json
from io import StringIO
import pandas as pd
from itertools import chain
import os
from django.conf import settings
import clusterCAD.celery as celery

def search(request):
    timeTaken = 0
    if request.method != 'POST':
        # if not a post, show the search form

        if 'aainput' in request.GET: 
            # prefill the search form with a referring AA sequence
            aainput = urlunquote(request.GET['aainput'])
            messages.success(request, 'Imported AA sequence from previous page')
        elif 'subunit' in request.GET: 
            # prefill the search form with a referring PKS subunit id
            subunitid = int(request.GET['subunit'])
            assert 0 <= subunitid
            subunit = Subunit.objects.get(id=subunitid)
            aainput = subunit.sequence
            messages.success(request, 
                'Imported sequence for PKS cluster %s subunit %s' 
                % (subunit.cluster.description, subunit.name))
        else:
            aainput = ''
        context = {'aainput': aainput}
        return render(request, 'sequencesearch.html', context)

    elif 'aainput' in request.POST:
        # if this is a post, validate the input and
        # run alignment

        # validate all inputs
        searchDatabase = str(request.POST['searchDatabase'])
        maxHits = int(request.POST['maxHits'])
        assert 1 <= maxHits <= 10000
        inputs = request.POST['aainput']
        evalue = float(request.POST['evalue'])
        assert 0.0 <= evalue <= 10.0
        
        showAllDomains = int(request.POST['showAllDomains'])
        assert 0 <= showAllDomains <= 1
        inputs = inputs.strip()
        if len(re.findall('>.*?\n', inputs)) >= 2:
            message = 'Error: Multiple queries detected, please remove until only one query is present'
            return JsonResponse({"errormessage": message}, status=406)
        inputs = re.sub('^>.*?\n', '', inputs)
        inputs = re.sub('\s', '', inputs)
        if len(inputs) == 0:
            message = 'Error: Empty Query'
            return JsonResponse({"errormessage": message}, status=406)
        if len(inputs) > 50000:
            message = 'Error: max query size is 50,000 residues'
            return JsonResponse({"errormessage": message}, status=406)

        # run Blast 
        if searchDatabase == "reviewed":
            db=os.path.join(
                    settings.BASE_DIR, 
                    'pipeline', 'data', 'blast', 'clustercad_subunits_reviewed',
                )
        else:
            db=os.path.join(
                    settings.BASE_DIR, 
                    'pipeline', 'data', 'blast', 'clustercad_subunits_all',
                )

        context = {
            'queryResidues': len(inputs),
            'evalue': str(evalue),
            'maxHits': str(maxHits),
            'showAllDomains': showAllDomains,
            'searchDatabase': str(searchDatabase),
        }

        # launch search as a celery task
        resultTask = blast.delay(query=inputs, 
            evalue=evalue, max_target_seqs=maxHits, sortOutput=True, database=db, context=context)

        task_id = resultTask.id

        # send the running task id back to the frontend
        return JsonResponse({"task_id": task_id}, status=202)
    
    else:
        raise Http404

def status(request, taskid):
    result = blast.AsyncResult(taskid)
    return JsonResponse({"task_status": result.state}, status=200)

def results(request, taskid):

    result = blast.AsyncResult(taskid)

    try:
        results, timeTaken, queries_found, submitcontext = result.get(timeout=30)
    except TimeoutError:
        messages.error(request, 
                'Job missing or timed out. Please try again. Reported status: ' + result.state)
        return render(request, 'sequencesearch.html')

    showAllDomains = submitcontext['showAllDomains']

    #If no queries found, then no hits
    if not queries_found:
        messages.error(request, 'No hits - please refine query!')
        return render(request, 'sequencesearch.html')

    #Read from jsonized pandas dataframe
    df = pd.read_json(results)
    
    #Get the subunits and process domains and modules
    subjectacc_separater = lambda x: Subunit.objects.get(id=x.split('_')[1], cluster__mibigAccession=x.split('_')[0])
    df['subunit'] = df['subject acc.ver'].map(subjectacc_separater)

    domain_getter = lambda x: Domain.objects.filter(module__subunit=x['subunit'], stop__gte=x['s. start'], start__lte=x['s. end']).select_subclasses().order_by('start')
    df['domains'] = df.apply(domain_getter, axis=1)

    module_getter = lambda x: list(set([domain.module for domain in x['domains']]))
    df['modules'] = df.apply(module_getter, axis=1)

    
    #User option to show all domains
    if showAllDomains:
        show_all_domains_func = lambda x: [{'module': module, 'domains': list(Domain.objects.filter(module=module).select_subclasses().order_by('start'))} for module in x['modules']]
        df['modules'] = df.apply(show_all_domains_func, axis=1)

    #Other option to show only relevant domains
    else:
        module_dict_getter = lambda x: [{'module': module, 'domains': list(x['domains'].filter(module=module))} for module in x['modules']]
        df['modules'] = df.apply(module_dict_getter, axis=1)

    #Keep useful columns only
    df = df[['subunit','modules','q. start','q. end','s. start','s. end','evalue','bit score','domains']]
    
    df = df.sort_values('bit score', ascending=False)
    alignments = df

    #Set cache
    # cache.set((inputs, evalue, maxHits, showAllDomains, searchDatabase), alignments,
    #     60 * 60 * 24 * 7) # cache for one week

    # get domain options to display in UI
    get_all_domains = lambda x: [domain for module in x['modules'] for domain in module['domains']]
    domains = alignments.apply(get_all_domains, axis=1).tolist()

    #flatten the list
    domains = list(chain.from_iterable(domains))

    #Extract desired types for html showing
    ats = list(filter(lambda d: isinstance(d, AT), domains))
    atsubstrates = list(set([at.substrate for at in ats]))
    krs = list(filter(lambda d: isinstance(d, KR), domains))
    krtypes = list(set([kr.type for kr in krs]))
    cs = list(filter(lambda d: isinstance(d, C), domains))
    ctypes = list(set([c.type for c in cs]))
    rs = list(filter(lambda d: isinstance(d, R), domains))
    rtypes = list(set([r.type for r in rs]))
    boolDomains = []
    for domain in (DH, ER, cMT, oMT, E, Cy, nMT, F, AOX, X):
        thesedomains = list(filter(lambda d: isinstance(d, domain), domains))
        typelist = list(set([str(d) for d in thesedomains]))
        if len(typelist) > 0:
            boolDomains.append((domain.__name__, typelist))
    tes = list(filter(lambda d: isinstance(d, TE), domains))
    tetypes = list(set([str(te) for te in tes]))
    adomains = list(filter(lambda d: isinstance(d, A), domains))
    asubstrates = list(set([a.substrate for a in adomains]))
    alignments = alignments.values.tolist()

    # create context dict of all results to show the user
    context = {
        'alignments': alignments,
        'timeTaken': timeTaken,
        'atsubstrates': atsubstrates,
        'asubstrates': asubstrates,
        'krtypes': krtypes,
        'ctypes': ctypes,
        'rtypes': rtypes,
        'boolDomains': boolDomains,
        'tetypes': tetypes,
    }

    # merge in submit context entries
    context.update(submitcontext)
    
    return render(request, 'sequenceresult.html', context)    
    
