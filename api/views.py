from django.shortcuts import render
from django.http import HttpResponse
from django.http import Http404
from json import dumps

def api(request):
    # if request.is_ajax():
        input = request.GET['integer']
        
        response = {
            'newinteger': int(input) + 1,
        }
        
        return HttpResponse(dumps(response), 'text/json')

    # else:
    #     # raise Http404
    #     return HttpResponse('Error: not an AJAX request!', 'text/json')
