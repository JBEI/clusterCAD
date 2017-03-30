import re
from django.shortcuts import render
from model_utils.managers import InheritanceManager
from django.http import HttpResponse
from .models import Cluster, Module, Subunit, Domain
from django.http import Http404
from json import dumps

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

def domainLookup(request):
    if request.is_ajax():
        try:
            domainid = request.GET['domainid']
            domain = Domain.objects.filter(id=int(domainid)).select_subclasses()[0]
            response = {
                'name': '%s subunit %s module %s: %s domain' % (domain.module.subunit.cluster.description, domain.module.subunit.name, domain.module.order, repr(domain)),
                'start': str(domain.start),
                'stop': str(domain.stop),
                'annotations': str(domain),
                'AAsequence': domain.getAminoAcidSequence(),
                'DNAsequence': domain.getNucleotideSequence(),
            }
        except:
            raise Http404
        return HttpResponse(dumps(response), 'text/json')
    else:
        raise Http404

def subunitLookup(request):
    if request.is_ajax():
        try:
            subunitid = request.GET['subunitid']
            subunit = Subunit.objects.get(id=int(subunitid))
            response = {
                'name': '%s subunit %s' % (subunit.cluster.description, subunit.name),
                'start': str(subunit.start),
                'stop': str(subunit.stop),
                'genbankAccession': subunit.genbankAccession,
                'genbankAccessionShort': re.sub("\.\d+$", "", subunit.genbankAccession),
                'AAsequence': subunit.getAminoAcidSequence(),
                'DNAsequence': subunit.getNucleotideSequence(),
            }
        except:
            raise Http404
        return HttpResponse(dumps(response), 'text/json')
    else:
        raise Http404
