#!/usr/bin/python3

import os
import sys
sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models
import domainLevelSearch.models

def isActive(domain):
    # returns the catalytic activity of a domain as a boolean

    if type(domain) in {
                            pks.models.CAL,
                            pks.models.AT,
                            pks.models.TE,
                            pks.models.KS,
                            pks.models.ACP,
                            pks.models.PCP,
                       }:
        return True # these domains are always active but don't contain the boolean internally
    else:
        return domain.active

def moduleToTuple(thisModule, inactives=True):
    # export a PKS module as a list of text strings
    if inactives:
        domains = thisModule.domains()
    else:
        domains = filter(lambda d: isActive(d), thisModule.domains())
    return tuple(repr(domain) + '_' + str(domain) for domain in domains)

nrpClusters = {domain.module.subunit.cluster.mibigAccession for domain in pks.models.A.objects.all()}

for thispks in pks.models.Cluster.objects.exclude(mibigAccession__in=nrpClusters).filter(reviewed=True):

    print('creating architecture string for: ' + str(thispks))
    # for this PKS get a list of domain text strings
    pksarchitecture = [domainText
        for subunit in thispks.subunits()
        for module in subunit.modules()
        for domainText in moduleToTuple(module, inactives=False)]

    # convert domain text strings into single unicode chars, and collapse to string
    domainChars = [chr(domainLevelSearch.models.DomainChar.objects.get_or_create(domainString=domainString)[0].id)
        for domainString in pksarchitecture]
    archString = ''.join(domainChars)

    # store architecture string in the database
    ClusterStringObj, created = domainLevelSearch.models.ClusterString.objects.get_or_create(cluster = thispks)
    ClusterStringObj.archString = archString
    ClusterStringObj.save()

