from django.db import models
from rdkit import Chem as chem

class Compound(models.Model):
    inchiKey = models.CharField(max_length=27, primary_key=True)
    name = models.CharField(max_length=256)
    smiles = models.TextField()

    def computeInchi(self):
        self.inchiKey = chem.InchiToInchiKey(chem.MolToInchi(chem.MolFromSmiles(self.smiles)))  
        return self.inchiKey
