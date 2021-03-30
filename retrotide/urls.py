# Based on example from medium article
# https://medium.com/codex/deploying-react-through-djangos-static-files-part-1-dev-setup-8a3a7b93c809

from django.contrib import admin
from django.urls import path, re_path
from django.shortcuts import render

def render_react(request):
    return render(request, "retrotide_index.html")
  
urlpatterns = [
    re_path(r"^$", render_react),
    re_path(r"^(?:.*)/?$", render_react),
]
