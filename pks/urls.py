from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, {'show': 'reviewed'}),
    url(r'^all/?$', views.index, {'show': 'all'}),
    url(r'^domainLookup/?$', views.domainLookup, name='domainLookup'),
    url(r'^subunitLookup/?$', views.subunitLookup, name='subunitLookup'),
    url(r'^(?P<mibigAccession>\S+)/$', views.details, name='details'),
]
