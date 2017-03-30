from django.shortcuts import render
from django.http import HttpResponse
from .models import Cluster, Module, Subunit
from django.http import Http404

def index(request):
    try:
        clusters=Cluster.objects.order_by('description')
    except Cluster.DoesNotExist:
        raise Http404

    clusterlist = [] 
    for cluster in clusters:
        clusterDict = {
            'clusterObject': cluster,
            'subunitCount': Subunit.objects.filter(cluster=cluster).count(),
            'moduleCount': Module.objects.filter(subunit__cluster=cluster).count(),
        }
        clusterlist.append(clusterDict)

    context={'clusters': clusterlist}

    return render(request, 'index.html', context)

def details(request, mibigAccession):
    try:
        cluster=Cluster.objects.get(mibigAccession=mibigAccession)
    except Cluster.DoesNotExist:
        raise Http404

    if 'mark' in request.GET:
        mark = int(request.GET['mark'])
    else:
        mark = -1 

    context={
            'cluster': cluster, 
            'architecture': cluster.architecture(),
            'mark': mark,
            'notips': ('KS', 'ACP', 'PCP'),
    }

    return render(request, 'details.html', context)
