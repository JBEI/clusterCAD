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
def moduleNumber(module):
    pos = list(Module.objects.filter(subunit__cluster = module.subunit.cluster).order_by('subunit__order', 'order')).index(module)
    return str(pos + 1)

@register.filter
def urlq(str):
    return urlquote(str)
