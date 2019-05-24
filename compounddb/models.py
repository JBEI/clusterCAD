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
        # returns the Inchi key for this compound
        self.inchiKey = chem.InchiToInchiKey(chem.MolToInchi(chem.MolFromSmiles(self.smiles)))  
        return self.inchiKey

    def mol(self):
        # returns and RDKit mol object for this compound
        return chem.MolFromSmiles(self.smiles)

    def aMW(self):
        # returns the Average Molecular Weight
        with connection.cursor() as c:
            c.execute('SELECT mol_amw(m) FROM compounddb_compound '\
                      'WHERE "inchiKey" = \'%s\';' \
                      % self.inchiKey)
            return c.fetchone()[0]

    def save(self, *args, **kwargs):
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

        # add fingerprints to fingerprint table
        addFingerprint(self.inchiKey)

    @classmethod
    def atomPairSearch(cls, querySmiles, maxHits=10, minSim=0.5, reviewedOnly=False):
        # perform atom pair similarity search
        # result is a list of tuples where each tuple is of the form
        # (tanimito similarity as integer, Compound object)

        with connection.cursor() as c:
            c.execute('SET rdkit.tanimoto_threshold=%s;' % str(minSim))
            if reviewedOnly:
                # get only compounds which we can join to a reviewed PKS gene cluster
                c.execute(
                    'SELECT DISTINCT ON (similarity, "inchiKey") '
                    'similarity, "inchiKey" '
                    'FROM get_ap_neighbors(\'%s\'), pks_module, pks_subunit, pks_cluster '
                    'WHERE "inchiKey"=pks_module.product_id AND '
                    'pks_module.subunit_id=pks_subunit.id AND '
                    'pks_subunit.cluster_id=pks_cluster."mibigAccession" AND '
                    'pks_cluster.reviewed=TRUE '
                    'ORDER BY similarity DESC LIMIT %s;'
                    % (querySmiles, str(maxHits))
                )
                    
            else:
                # get all compounds
                c.execute('SELECT similarity, "inchiKey" FROM get_ap_neighbors(\'%s\') LIMIT %s;' % (querySmiles, str(maxHits)))
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

def addFingerprint(inchiKey):
    # add fingerprint of a compound with the given incihKey to the fingerprint table.
    # Compound must be already in the compounddb_compound table.
    with connection.cursor() as c:
        c.execute('SELECT "inchiKey" FROM rdk.fps WHERE "inchiKey" = \'%s\' LIMIT 1;' % inchiKey)
        if c.fetchone():
            return False
        c.execute('INSERT INTO rdk.fps SELECT "inchiKey", atompair_fp(m) AS apfp '\
                  'FROM compounddb_compound '\
                  'WHERE "inchiKey" = \'%s\'' % inchiKey)
