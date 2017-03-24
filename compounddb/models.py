from django.db import models
from rdkit import Chem as chem
from django.db import connection
from django.db.models.signals import pre_delete
from django.dispatch import receiver

class Compound(models.Model):
    # this class stores a small molecule structure

    inchiKey = models.CharField(max_length=27, primary_key=True)
    name = models.CharField(max_length=256)
    smiles = models.TextField()

    def computeInchiKey(self):
        self.inchiKey = chem.InchiToInchiKey(chem.MolToInchi(chem.MolFromSmiles(self.smiles)))  
        return self.inchiKey

    def mol(self):
        # returns and RDKit mol object for this compound
        return chem.MolFromSmiles(self.smiles)

    def save(self, *args, **kwargs):
        # if this is the first save, add mols column to database
        addMolsCol()

        # compute inchiKey from structure
        if not self.inchiKey:
            self.computeInchiKey()
        super(Compound, self).save(*args, **kwargs)

        # add mol object to database column
        with connection.cursor() as c:
            c.execute('UPDATE compounddb_compound '\
                      'SET m = mol_from_smiles(\'%s\') '\
                      'WHERE "inchiKey" = \'%s\';' \
                      % (self.smiles, self.inchiKey))

        # add fingerprints
        addFingerprint(self.inchiKey)

    @classmethod
    def atomPairSearch(cls, querySmiles, maxHits=10, minSim=0.5):
        # perform atom pair similarity search
        # result is a list of tuples where each tuple is of the form
        # (tanimito similarity as integer, Compound object)

        # validate query
        try:
            query = chem.MolFromSmiles(querySmiles)
            chem.SanitizeMol(query)
        except:
            raise ValueError('Invalid input smiles') 

        with connection.cursor() as c:
            c.execute('SET rdkit.tanimoto_threshold=%s;' % str(minSim))
            c.execute('SELECT similarity, "inchiKey" FROM get_ap_neighbors(\'%s\') limit %s;' % (querySmiles, str(maxHits)))
            results = c.fetchall()
        objectResult = []
        for result in results:
            objectResult.append((result[0], cls.objects.get(inchiKey = result[1])))
        return objectResult

@receiver(pre_delete, sender=Compound)
def delete_fp(sender, instance, **kwargs):
    # drop compound from fingerprint table if deleted
    with connection.cursor() as c:
        c.execute('DELETE FROM rdk.fps WHERE "inchiKey" = \'%s\';' % instance.inchiKey)

def addMolsCol():
    # add RDKit mol column to database table

    # if column already exists, return false
    with connection.cursor() as cursor:
        cursor.execute('SELECT * FROM compounddb_compound LIMIT 0;')
        colnames = [desc[0] for desc in cursor.description]
    if 'm' in colnames:
        return False

    # add indexed mol column
    with connection.cursor() as c:
        c.execute('CREATE EXTENSION IF NOT EXISTS rdkit;')
        c.execute('CREATE SCHEMA IF NOT EXISTS rdk;')
        c.execute('ALTER TABLE compounddb_compound '\
                       'ADD COLUMN m mol;') 
        c.execute('CREATE INDEX molidx ON compounddb_compound USING gist (m);')
    return True

def createFingerprintTable():
    # add atom pair fingerprint table to database

    with connection.cursor() as c:
        c.execute('CREATE TABLE IF NOT EXISTS rdk.fps '\
            '("inchiKey" character varying(27) NOT NULL, '\
            'apfp sfp);')
        c.execute('ALTER TABLE ONLY rdk.fps '\
            'ADD PRIMARY KEY ("inchiKey");')
        c.execute('CREATE INDEX fps_apfp_idx ON rdk.fps USING gist (apfp);')
        
        # define atom pair search function
        c.execute('CREATE OR REPLACE FUNCTION get_ap_neighbors(smiles text) '\
                  'RETURNS TABLE("inchiKey" character, m mol, similarity double precision) AS '\
                  '$$ '\
                  'SELECT "inchiKey",m,tanimoto_sml(atompair_fp(mol_from_smiles($1::cstring)),apfp) AS similarity '\
                  'FROM rdk.fps JOIN compounddb_compound USING ("inchiKey") '\
                  'WHERE atompair_fp(mol_from_smiles($1::cstring))%apfp '\
                  'ORDER BY similarity DESC; '\
                  '$$ LANGUAGE SQL STABLE;')

def addFingerprint(inchiKey):
    with connection.cursor() as c:
        c.execute('SELECT * FROM information_schema.tables WHERE table_schema = \'rdk\';')
        if not c.fetchone():
            createFingerprintTable()
        c.execute('SELECT "inchiKey" FROM rdk.fps WHERE "inchiKey" = \'%s\' LIMIT 1;' % inchiKey)
        if c.fetchone():
            return False
        c.execute('INSERT INTO rdk.fps SELECT "inchiKey", atompair_fp(m) AS apfp '\
                  'FROM compounddb_compound '\
                  'WHERE "inchiKey" = \'%s\'' % inchiKey)
