from django import template
from django.utils.http import urlquote
from rdkit import Chem
from pks.models import ACP, cMT, PCP, Module
import re

register = template.Library()

@register.filter
def classname(obj):
    return obj.__class__.__name__

@register.filter
def smiles(mol):
    if mol:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        return urlquote(smiles)
    else:
        return False 

@register.filter
def unquotedSmiles(mol):
    return Chem.MolToSmiles(mol, isomericSmiles=True)

@register.filter
def stripTrailingVersion(accession):
    return re.sub("\.\d+$", "", accession)

@register.filter
def acpDisplace(module):
    displace = 10 + (module.count() * 35)
    typeList = [type(x) for x in module]
    if not ACP in typeList:
        displace += 40
    if cMT in typeList:
        displace += 10 
    if PCP in typeList:
        displace += 10 
    return displace

@register.filter
def moduleNumber(module):
    pos = list(Module.objects.filter(subunit__cluster = module.subunit.cluster).order_by('subunit__order', 'order')).index(module)
    return str(pos + 1)
