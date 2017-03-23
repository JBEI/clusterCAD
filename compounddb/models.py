from django.db import models

class Compound(models.Model):
    inchiKey = models.CharField(max_length=27, primary_key=True)
    name = models.CharField(max_length=256)
    smiles = models.TextField()
