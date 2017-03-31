from django import template
from django.template.defaultfilters import stringfilter
from django.utils.http import urlquote

register = template.Library()

@register.filter
def urlq(str):
    return urlquote(str)
