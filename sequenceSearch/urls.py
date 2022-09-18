from django.conf.urls import url
from django.urls import path

from . import views

urlpatterns = [
    path('status/<slug:taskid>/', views.status),
    path('results/<slug:taskid>/', views.results),
    url(r'^$', views.search),
]
