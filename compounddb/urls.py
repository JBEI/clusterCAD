from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^renderer/(\S+)/(\d{1,4})$', views.renderer),
    url(r'^mcsrenderernew/(?P<smiles1>\S+)_(?P<smiles2>\S+)/(?P<width>\d{1,4})/(?P<align>\d{1})$', views.mcsrenderer),
    url(r'^mcsrenderer/(\S+)_(\S+)_(\S+)/(\d{1,4})$', views.mcsrenderer),
]
