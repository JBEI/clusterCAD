from django.shortcuts import render
from django.http import HttpResponse
from .models import Cluster
from django.http import Http404

def index(request):
    try:
        clusters=Cluster.objects.order_by('description')
    except Cluster.DoesNotExist:
        raise Http404

    context={'clusters': clusters}

    return render(request, 'index.html', context)

def details(request, pksId):
    try:
        cluster=Cluster.objects.get(id__iexact=pksId)
    except Cluster.DoesNotExist:
        raise Http404

    context={'cluster': cluster, 'architecture': cluster.architecture(products=True)}

    return render(request, 'details.html', context)
