from django.db import models
from rdkit import Chem as chem
from django.db import connection

class Compound(models.Model):
    inchiKey = models.CharField(max_length=27, primary_key=True)
    name = models.CharField(max_length=256)
    smiles = models.TextField()

    def computeInchi(self):
        self.inchiKey = chem.InchiToInchiKey(chem.MolToInchi(chem.MolFromSmiles(self.smiles)))  
        return self.inchiKey

    def save(self, *args, **kwargs):
        createMolsTable()
        if not self.inchiKey:
            self.computeInchi()
        with connection.cursor() as cursor:
            cursor.execute('DELETE FROM mols WHERE "inchiKey" = \'%s\';' % self.inchiKey)
            cursor.execute('INSERT INTO mols ("inchiKey", m) '\
                           'VALUES (\'%s\', mol_from_smiles(\'%s\'));' % (self.inchiKey, self.smiles))
        super(Compound, self).save(*args, **kwargs)

def createMolsTable():
    if 'mols' in connection.introspection.table_names():
        return False
    with connection.cursor() as cursor:
        cursor.execute('CREATE EXTENSION IF NOT EXISTS rdkit;')
        cursor.execute('CREATE SCHEMA IF NOT EXISTS rdk;')
        cursor.execute('CREATE TABLE IF NOT EXISTS mols '\
            '("inchiKey" character varying(27) NOT NULL, '\
            'm mol '\
            ');')
        cursor.execute('ALTER TABLE ONLY mols '\
            'ADD CONSTRAINT mols_pkey PRIMARY KEY ("inchiKey");')
        cursor.execute('CREATE INDEX molidx ON mols USING gist (m);')
    return True
