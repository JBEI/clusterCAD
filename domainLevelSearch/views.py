from django.shortcuts import render
from django.http import HttpResponse
from django.http import Http404
from json import loads, dumps
import logging

def domainSearch(request):

    # parse input JSON
    input = loads(request.GET['modules[]'])
    domainList = input['domainList']

    # convert input JSON into a parseable python object

    return HttpResponse(dumps(domainList), 'text')

    # context = {}
    
    # return render(request, 'domainsearchresults.html', context)

    # else:
    #     # raise Http404
    #     return HttpResponse('Error: not an AJAX request!', 'text/json')
