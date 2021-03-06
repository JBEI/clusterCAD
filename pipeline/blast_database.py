#!/usr/bin/python3

import os, sys
from Bio import Seq, SeqRecord, SeqIO
from Bio.Alphabet import IUPAC
import argparse

def main():
  parser = argparse.ArgumentParser(description='blast database creation')
  parser.add_argument('--name', type=str, default='clustercad_subunits_reviewed', metavar='N',
                        help='name for clustercad database')
  parser.add_argument('--include_all', type=int, default=0, metavar='N',
                        help='0 for reviewed only, 1 for all')
  args = parser.parse_args()


  sys.path.insert(0, '/clusterCAD')
  os.environ.setdefault("DJANGO_SETTINGS_MODULE", "clusterCAD.settings")
  import django
  django.setup()
  import pks.models
  from django.conf import settings

  # resolve the path to the blast database (reviewed)
  dbPath = os.path.join(
      settings.BASE_DIR, 
      'pipeline', 'data', 'blast', args.name,
  )



  # loop over all AA sequences in ClusterCAD and write them to the database
  # file
  with open(dbPath, 'w') as f:
      if args.include_all == 0:
        for cluster in pks.models.Cluster.objects.filter(reviewed=True):
            for subunit in cluster.subunits():
                modacc = cluster.mibigAccession + '_' + str(subunit.id) 
                sseq = SeqRecord.SeqRecord(Seq.Seq(subunit.sequence, IUPAC.protein),
                                           id=modacc,
                                           name=subunit.name,
                                           description=''
                                          )
                SeqIO.write(sseq, f, "fasta")
      else:
        for cluster in pks.models.Cluster.objects.all():
            for subunit in cluster.subunits():
                modacc = cluster.mibigAccession + '_' + str(subunit.id) 
                sseq = SeqRecord.SeqRecord(Seq.Seq(subunit.sequence, IUPAC.protein),
                                           id=modacc,
                                           name=subunit.name,
                                           description=''
                                          )
                SeqIO.write(sseq, f, "fasta")


if __name__ == '__main__':
    main()
