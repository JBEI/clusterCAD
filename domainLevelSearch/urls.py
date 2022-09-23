from django.conf.urls import url

from . import views
from clusterCAD.urls import CsrfView

urlpatterns = [
    url(r'^query/$', CsrfView.as_view(template_name='domainsearchquery.html')),
    url(r'^$', views.domainSearch),
]
