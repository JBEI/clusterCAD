# -*- coding: utf-8 -*-
# This migration sets up the RDKit postgresql catridge
# it was written by hand, and cannot be generated with makemigrations

from __future__ import unicode_literals
from django.db import migrations, models

class Migration(migrations.Migration):

    dependencies = [
        ('compounddb', '0001_initial'),
    ]

    operations = [
        migrations.RunSQL('CREATE EXTENSION IF NOT EXISTS rdkit;'),
        migrations.RunSQL('CREATE SCHEMA IF NOT EXISTS rdk;'),
        migrations.RunSQL('ALTER TABLE compounddb_compound '\
                          'ADD COLUMN m mol;'),
        migrations.RunSQL('CREATE INDEX molidx ON compounddb_compound USING gist (m);'),
        migrations.RunSQL('CREATE TABLE IF NOT EXISTS rdk.fps '\
                          '("inchiKey" character varying(27) NOT NULL, '\
                          'apfp sfp);'),
        migrations.RunSQL('ALTER TABLE ONLY rdk.fps '\
                          'ADD PRIMARY KEY ("inchiKey");'),
        migrations.RunSQL('CREATE INDEX fps_apfp_idx ON rdk.fps USING gist (apfp);'),
        migrations.RunSQL('CREATE OR REPLACE FUNCTION get_ap_neighbors(smiles text) '\
                          'RETURNS TABLE("inchiKey" character, similarity double precision) AS '\
                          '$$ '\
                          'SELECT "inchiKey",tanimoto_sml(atompair_fp(mol_from_smiles($1::cstring)),apfp) AS similarity '\
                          'FROM rdk.fps '\
                          'WHERE atompair_fp(mol_from_smiles($1::cstring))%apfp '\
                          'ORDER BY similarity DESC; '\
                          '$$ LANGUAGE SQL STABLE;'),
    ]
