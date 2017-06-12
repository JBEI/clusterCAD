#!/usr/bin/python3

import os, sys
from Bio import Seq, SeqRecord, SeqIO
from Bio.Alphabet import IUPAC

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

with open('/clusterCAD/pipeline/data/blast/clustercad_subunits', 'w') as f:
    for cluster in pks.models.Cluster.objects.all():
        for subunit in cluster.subunits():
            modacc = cluster.mibigAccession + '_' + str(subunit.id) 
            sseq = SeqRecord.SeqRecord(Seq.Seq(subunit.sequence, IUPAC.protein),
                                       id=modacc,
                                       name=subunit.name,
                                       description=''
                                      )
            SeqIO.write(sseq, f, "fasta")
