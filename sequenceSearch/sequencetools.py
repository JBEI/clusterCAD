import os
from io import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from django.conf import settings
from pks.models import Subunit, Domain
from model_utils.managers import InheritanceManager
from copy import deepcopy
from time import time
import json
from pks.models import AT, KR, DH, ER, cMT, oMT, TE, Subunit, Domain
import pandas as pd

def blast(
        query, 
        context,
        evalue=10.0, 
        max_target_seqs=10, 
        sortOutput=True,
        database=os.path.join(
                settings.BASE_DIR, 
                'pipeline', 'data', 'blast', 'clustercad_subunits_reviewed',
            ),

    ):
    # run blast and return results as a list
    # of alignment dicts with the following structure:
    # {'alignment': alignment, 'subunit': subunit, 'hsps': hsps}
    # where hsps have the following structure:
    # {'hsp': hsp, 'modules': modules}
    # and modules has the following structure:
    # {'module': module, 'domains': list(domains)} 

    # Check inputs
    assert isinstance(evalue, float)
    assert 0.0 <= evalue <= 10.0
    assert isinstance(max_target_seqs, int)
    assert 1 <= max_target_seqs <= 10000
    assert isinstance(query, str)

    # Convert query to fasta format 
    queryStringIO = StringIO()
    queryRecord = SeqRecord(Seq(query, IUPAC.protein), id='query')
    SeqIO.write(queryRecord, queryStringIO, "fasta")
    queryFasta = queryStringIO.getvalue()
    queryStringIO.close()

    # Column names returned by NcbiblastpCommandline:
    # query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

    # Run blast
    start = time()
    blastp_cline = NcbiblastpCommandline(
                                         db=database,
                                         evalue=evalue,
                                         outfmt=7,
                                         num_threads=2
                                         )
    result, stderr = blastp_cline(stdin=queryFasta)

    # Get column names from blast output
    results = result.split("\n")
    columns = results[3].split(":")[-1].split(",")
    columns = [i.lstrip() for i in columns]

    # If no hits found, return immediately
    if results[3] == "# 0 hits found":
        return None, None, False

    # Parse blast output and delete files
    # Remeber bitscore cant convert to int by rule 'safe'
    # Load results into a pandas dataframe
    resultIO = StringIO(result)
    df = pd.read_csv(resultIO, sep="\t", comment="#", names=columns)
    resultIO.close()


    # If sort outputs, then sort according to bitscore
    if sortOutput:
        df = df.sort_values('bit score', ascending=False)
        

    # Return only up to max_target_seq results from query
    if max_target_seqs < len(df):
        df = df[:max_target_seqs]

    end = time()

    return df.to_json(), str(int(end-start)), True, context

    
