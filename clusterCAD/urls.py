"""clusterCAD URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.10/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url, include
from django.contrib import admin
from django.views.generic import TemplateView
from django.views.decorators.csrf import ensure_csrf_cookie
from django.utils.decorators import method_decorator

# define a new csrf view for static pages to set the csrf token
# even on static pages
class CsrfView(TemplateView):
    
    @method_decorator(ensure_csrf_cookie)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^pks/', include('pks.urls')),
    url(r'^compounddb/', include('compounddb.urls')),
    url(r'^structureSearch/', include('structureSearch.urls')),
    url(r'^sequenceSearch/', include('sequenceSearch.urls')),
    url(r'^api/', include('domainLevelSearch.urls')),
    url(r'^$', CsrfView.as_view(template_name='home.html')),
    url(r'^about/$', CsrfView.as_view(template_name='about.html')),
] 
