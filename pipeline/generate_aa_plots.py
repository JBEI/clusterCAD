#!/usr/bin/python3

import matplotlib
matplotlib.use('Agg')

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import hashlib
import json
from io import StringIO

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
from django.core.files import File
from django.core.files.base import ContentFile

import pks.models

plotCache = 'data/aa_sequence_analysis/plotCache/'

def plot_heatmap(values, labels, mapping=None):
    """Generates a 1D heat map in lines of length n.
    
    Args:
      values: A list or string of values
      labels: A list or string of labels
    """

    # This is the number of positions to plot per line
    # Note: Spacing between lines will break for n <= 25 
    #       Sizing of the last line will break for n >= 100
    n = 50

    # Number of lines that will be added to plot
    l = int(np.ceil(len(labels) / float(n)))

    if mapping:
        l += 1

    # These multiplicative factors are important for keeping
    # the scaling of the last axis consistent with the others
    # scaling seems to get messed up when axis feels squished
    fig = plt.figure(figsize=(0.25 * n, 0.6 * l + 0.5), dpi=80)

    # For each line
    for i in range(l):
        if mapping and i == (l - 1):
            vs = list(mapping.values())
            ls = ''.join(mapping.keys())
        else:
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
        if mapping and i == (l - 1):
            ax.set_yticklabels(['KEY'], size='small')
        else:
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

        if i == l - 1 and not mapping:
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
    return outputSVG.getvalue()

plotless_subunits = pks.models.Subunit.objects.filter(accPlotFile='', ssPlotFile='')

# don't try to make plots for subunits missing acc20 or ss8 data
# also fill in missing acc20 or ss8 data as [DATA IS UNAVAILABLE]
# also don't try to make plots for those with acc20 or ss8 stating [DATA IS UNAVAILABLE]
subunits = []
for subunit in plotless_subunits:
    if subunit.acc20 == "" or subunit.ss8 == "":
        subunit.acc20 = "[DATA UNAVAILABLE]"
        subunit.acc = "[DATA UNAVAILABLE]" # for display on the interface
        subunit.ss8 = "[DATA UNAVAILABLE]"
        subunit.ss = "[DATA UNAVAILABLE]" # for display on the interface
        subunit.save()
        continue
    elif subunit.acc20 == "[DATA UNAVAILABLE]" or subunit.ss8 == "[DATA UNAVAILABLE]":
        continue
    else:
        subunits.append(subunit)

# exit if everything has been plot already
if len(subunits) == 0:
    sys.exit("Done generating all plots.")

# process only one subunit
subunit = subunits[0]

mapping = {'C':50, 
           'G':60, 'H':70, 'I':80, 
           'E':30, 'B':10, 
           'T':100, 
           'S':0}

# generate plot 
print("Loading plots for " + subunit.cluster.description + " subunit " + subunit.name)

# check if acc20/ss8 data exists, if so process it
try:
    # generate solvent accessibility plot
    intList = [int(i) for i in subunit.acc20.split(',')]
    aaseq = list(subunit.getAminoAcidSequence())
    assert len(intList) == len(aaseq)
    hashstring = json.dumps(tuple(intList)) + json.dumps(tuple(aaseq))
    hashstring = hashlib.sha224(hashstring.encode()).hexdigest()
    filename = plotCache + hashstring + '.svg'
    if os.path.isfile(filename):
        print('\tfound cached acc20 plot')
        with open(filename, 'r') as svgFile: 
            accPlot = svgFile.read()
    else:
        print('\tno aac20 plot found, generating')
        accPlot = plot_heatmap(intList, aaseq)
        with open(filename, 'w') as svgFile:
            svgFile.write(accPlot)
    accplotname = str(subunit.id) + '.svg'
    subunit.accPlotFile.save(name=accplotname, content=ContentFile(accPlot), save=True)

    # generate secondary structure plot
    f = lambda x: mapping[x]
    ss_seq_nums = list(map(f, subunit.ss8))

    assert len(ss_seq_nums) == len(aaseq)
    hashstring = json.dumps(tuple(ss_seq_nums)) + json.dumps(tuple(aaseq))
    hashstring = hashlib.sha224(hashstring.encode()).hexdigest()
    filename = plotCache + hashstring + '.svg'
    if os.path.isfile(filename):
        print('\tfound cached ss8 plot')
        with open(filename, 'r') as svgFile: 
            ssPlot = svgFile.read()
    else:
        print('\tno ss8 plot found, generating')
        ssPlot = plot_heatmap(ss_seq_nums, aaseq, mapping)
        with open(filename, 'w') as svgFile:
            svgFile.write(ssPlot)
    ssplotname = str(subunit.id) + '.svg'
    subunit.ssPlotFile.save(name=ssplotname, content=ContentFile(ssPlot), save=True)

# handle ValueError that arises when split is attempted but acc20 data does not exist
except ValueError:
    # for some NRPS subunits acc20 data is not available. then ignore both acc20/ss8 and continue to next subunit
    print("\tacc20 (solvent accessibility) and/or ss8 (secondary structure) data not available. ")
    subunit.acc = "[DATA UNAVAILABLE]"
    subunit.ss = "[DATA UNAVAILABLE]"
    subunit.acc20 = "[DATA UNAVAILABLE]"
    subunit.ss8 = "[DATA UNAVAILABLE]"
    subunit.save()

# handle AssertionError that arises when sequence length is no longer equal to length of solvent accessibility values
except AssertionError:
    # may happen when mibig is updated with new sequences (rare). then ignore both acc20/ss8 and continue to next subunit
    print("\tAA sequence length does not match acc20 (solvent accessibility) and/or ss8 (secondary structure) value length.")
    subunit.acc = "[DATA UNAVAILABLE]"
    subunit.ss = "[DATA UNAVAILABLE]"
    subunit.acc20 = "[DATA UNAVAILABLE]"
    subunit.ss8 = "[DATA UNAVAILABLE]"
    subunit.save()
