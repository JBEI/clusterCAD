from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^renderer/(?P<smiles>\S*)/(?P<smiles2>\S*)$', views.renderer),
]
