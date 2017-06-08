import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
from django.shortcuts import render
from model_utils.managers import InheritanceManager
from django.http import HttpResponse
from .models import Cluster, Module, Subunit, Domain
from django.http import Http404
from json import dumps
from rdkit import Chem as chem
from io import StringIO
from django.core.cache import cache

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
        mark = [int(m) for m in request.GET['mark'].split(',')]
    else:
        mark = [] 

    architecture = cluster.architecture()

    # compute MCS percentages
    knownProduct = chem.MolFromSmiles(cluster.knownProductSmiles)
    mcs = chem.MolFromSmarts(cluster.knownProductMCS) 
    knownProductPercent = 100.0 * float(mcs.GetNumAtoms()) / float(knownProduct.GetNumAtoms())
    predictedProduct = architecture[-1][1][-1][0].product.mol()
    predictedProductPercent = 100.0 * float(mcs.GetNumAtoms()) / float(predictedProduct.GetNumAtoms())

    context={
            'cluster': cluster, 
            'architecture': architecture,
            'mark': mark,
            'notips': ('KS', 'ACP', 'PCP'),
            'predictedProductPercent': predictedProductPercent,
            'knownProductPercent': knownProductPercent,
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
                'id': str(subunit.id),
                'start': str(subunit.start),
                'stop': str(subunit.stop),
                'genbankAccession': subunit.genbankAccession,
                'genbankAccessionShort': re.sub("\.\d+$", "", subunit.genbankAccession),
                'AAsequence': subunit.getAminoAcidSequence(),
                'DNAsequence': subunit.getNucleotideSequence(),
                'ss': subunit.ss,
                'acc': subunit.acc, 
            }
        except:
            raise Http404
        return HttpResponse(dumps(response), 'text/json')
    else:
        raise Http404

def solventAccessibilityPlot(request, subunit):
    try:
        subunit=Subunit.objects.get(id=subunit)
    except Subunit.DoesNotExist:
        raise Http404
    intList = [int(i) for i in subunit.acc20.split(',')]
    aaseq = list(subunit.getAminoAcidSequence())
    assert len(intList) == len(aaseq)
    plot = plot_heatmap(intList, aaseq)
    return HttpResponse(plot, content_type='image/svg+xml')

def secondaryStruturePlot(request, subunit):
    try:
        subunit=Subunit.objects.get(id=subunit)
    except Subunit.DoesNotExist:
        raise Http404
    mapping = {'C':50, 
               'G':60, 'H':70, 'I':80, 
               'E':30, 'B':10, 
               'T':100, 
               'S':0}
    f = lambda x: mapping[x]
    ss_seq_nums = list(map(f, subunit.ss))
    aaseq = list(subunit.getAminoAcidSequence())
    assert len(aaseq) == len(aaseq)
    plot = plot_heatmap(ss_seq_nums, aaseq)
    return HttpResponse(plot, content_type='image/svg+xml')

def plot_heatmap(values, labels):
    """Generates a 1D heat map in lines of length n.
    
    Args:
      values: A list or string of values
      labels: A list or string of labels
    """
    cacheResults = cache.get(('plot_heatmap', values, labels))
    if cacheResults:
        return cacheResults

    # This is the number of positions to plot per line
    # Note: Spacing between lines will break for n <= 25 
    #       Sizing of the last line will break for n >= 100
    n = 50

    # Number of lines that will be added to plot
    l = int(np.ceil(len(labels) / float(n)))

    # These multiplicative factors are important for keeping
    # the scaling of the last axis consistent with the others
    # scaling seems to get messed up when axis feels squished
    fig = plt.figure(figsize=(0.25 * n, 0.6 * l + 0.5), dpi=80)

    # For each line
    for i in range(l):
        try:
            # Obtain values and labels for that line
            vs = values[i * n: (i + 1) * n]
            ls = labels[i * n: (i + 1) * n]
        except:
            vs = values[i * n: -1]
            ls = labels[i * n: -1]
        
        # Generate masked array from data
        mat = np.ma.array([vs, vs])

        # Add subplot for a single line
        ax = fig.add_subplot(l, 1, i + 1)
        cmap = plt.cm.RdBu
        # Deals with missing data appropriate
        # Specified using masked array (this is the purpose of the masked array)
        cmap.set_bad(color='0.5', alpha=1.0)
        pc = plt.pcolormesh(mat, cmap=cmap,
                            vmin=0, vmax=100)
        
        # Axis
        ax.set_xlim([0, n])
        ax.set_ylim([0, 1])
        plt.setp(ax.get_yticklines(), visible=False)
        ax.yaxis.tick_left()
        ax.set_yticks([0.5])
        ax.set_yticklabels([str(i * n + 1) + ' '], size='small')
        ax.xaxis.tick_bottom()
        ax.set_xticks(np.arange(0.5, len(ls) + 0.5, 1))
        ax.set_xticklabels(ls, size='x-small')
        ax.set_aspect('equal')
        # Set anchor to the left
        ax.set_anchor('W')
        
        # Tight layout must be declared before changing position of last axis
        plt.tight_layout()
        # Tight layout will also not leave space at the bottom for colorbar
        # Parameter set below can be used to fix spacing issues
        plt.subplots_adjust(left=0.05, right=0.95, bottom=0.20, top=0.95)

        if i == l - 1:
            # Must come after other settings to size last axis correctly
            ax.set_xlim([0, len(vs)], auto=False)
            # Not that get_position and set_position refer
            # to [left bottom width height] as fractions of plotting area
            pos = ax.get_position().get_points()
            # rescale to accommodate shortened sequence
            # ax.set_position([pos[0][0], pos[0][1],
            #                  pos[1][0], pos[1][1]])
            # add_axes refers to [left bottom width height]
            # as fractions of plotting area
            ax_c = fig.add_axes([pos[0][0], 0.1, 0.8, 0.03])
            plt.colorbar(pc, orientation="horizontal", cax=ax_c)

    outputSVG = StringIO()
    plt.savefig(outputSVG, format="svg")
    plotSVGstring = outputSVG.getvalue()

    # cache result forever
    cache.set(('plot_heatmap', values, labels), plotSVGstring, None)
    return plotSVGstring
