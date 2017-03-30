from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^domainLookup/?$', views.domainLookup, name='domainLookup'),
    url(r'^(?P<mibigAccession>\S+)/$', views.details, name='details'),
]
