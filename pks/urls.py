from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^solventAccessibilityPlot/(?P<subunit>\S+)/$', views.solventAccessibilityPlot, name='solventAccessibilityPlot'),
    url(r'^secondaryStruturePlot/(?P<subunit>\S+)/$', views.secondaryStruturePlot, name='secondaryStruturePlot'),
    url(r'^domainLookup/?$', views.domainLookup, name='domainLookup'),
    url(r'^subunitLookup/?$', views.subunitLookup, name='subunitLookup'),
    url(r'^(?P<mibigAccession>\S+)/$', views.details, name='details'),
]
