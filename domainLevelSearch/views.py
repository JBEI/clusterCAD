from django.shortcuts import render
from django.http import HttpResponse
from django.http import Http404
from json import loads, dumps
import Levenshtein
import domainLevelSearch.models
from time import time

def domainSearch(request):

    # parse JSON from POST data
    querydata = loads(request.body)

    # define module order for parsing POST data
    domainOrder = ['KS', 'AT', 'DH', 'ER', 'KR', 'ACP']
    moduleNo = 0
    queryList = []

    for module in querydata['params']['modules']:
        DomainListJson = module['DomainList']

        for thisdomain in domainOrder:

            # skip domains not listed in this module
            if thisdomain not in DomainListJson:
                continue

            # skip domains not present in this module
            if not DomainListJson[thisdomain]['present']:
                continue

            domainstring = thisdomain

            # add custom text for each domain type

            if thisdomain == 'KS':
                domainstring += '_domain'

            if thisdomain == 'AT':
                domainstring += '_substrate '
                domainstring += DomainListJson[thisdomain]['options']['selected']

                if moduleNo == 0:
                    domainstring += ', loading'
                else:
                    domainstring += ', non-loading'

            if thisdomain in {'DH', 'ER'}:
                domainstring += '_active'

            if thisdomain == 'KR':
                domainstring += '_type '
                domainstring += DomainListJson[thisdomain]['options']['selected']
                domainstring += ', active'

            if thisdomain == 'ACP':
                domainstring += '_substrate None'

            queryList.append(domainstring)

        moduleNo += 1

    # compile the list of domains into a query string

    queryString = ''
    for domainString in queryList:
        try:
            # see if it's in the database
            charobject = domainLevelSearch.models.DomainChar.objects.get(domainString=domainString)
        except DoesNotExist:
            # if this domain isn't in the database, set to 0
            charobject = chr(0)

        queryString += chr(charobject.id)

    # run search and compute time to run
    start = time()
    results = alignAll(queryString, tryTruncations=True)
    end = time()

    context = {
        'results': results,
        'timeTaken': "{:.2f}".format(end-start),
    }
    
    return render(request, 'domainsearchresults.html', context)

def alignAll(query, tryTruncations=False):
    # align query against all PKSs and return a sorted tuple of (mibig, editDistance)
    results = []
    for clusterstring in domainLevelSearch.models.ClusterString.objects.all():
        key = clusterstring.cluster
        value = clusterstring.archString
        bestResult = (key, Levenshtein.distance(query, value), 0)

        if tryTruncations:
            for truncation in range(1,len(value)):
                newResult = (key, Levenshtein.distance(query, value[0:-truncation]), truncation)
                if newResult[1] < bestResult[1]:
                    bestResult = newResult
        
        results.append(bestResult)
        
    return sorted(results, key=lambda x: x[1])

