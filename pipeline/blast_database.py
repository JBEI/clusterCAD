#!/usr/bin/python

import os, sys
from Bio import Seq, SeqRecord, SeqIO
from Bio.Alphabet import IUPAC

sys.path.insert(0, '/clusterCAD')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
import django
django.setup()
import pks.models

with open('./data/blast/subunits.fasta', 'wb') as f:
    for cluster in pks.models.Cluster.objects.all():
        index = 0
        for subunit in cluster.subunits():
            acc = subunit.genbankAccession.split('.')
            modacc = acc[0] + '.' + str(index + int(acc[1]))
            index += 1
            sseq = SeqRecord.SeqRecord(Seq.Seq(subunit.sequence, IUPAC.protein),
                                       id='gb|' + modacc,
                                       name=subunit.name,
                                       description=''
                                      )
            SeqIO.write(sseq, f, "fasta")
