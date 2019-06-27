import os
from io import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from pks.models import Subunit, Domain
from model_utils.managers import InheritanceManager
from copy import deepcopy
from time import time
from django.conf import settings

def blast(
        query, 
        evalue=10.0, 
        max_target_seqs=10, 
        sortOutput=True,
        database='reviewed'

    ):
    # run blast and return results as a list
    # of alignment dicts with the following structure:
    # {'alignment': alignment, 'subunit': subunit, 'hsps': hsps}
    # where hsps have the following structure:
    # {'hsp': hsp, 'modules': modules}
    # and modules has the following structure:
    # {'module': module, 'domains': list(domains)} 

    # check inputs
    assert isinstance(evalue, float)
    assert 0.0 <= evalue <= 10.0
    assert isinstance(max_target_seqs, int)
    assert 1 <= max_target_seqs <= 10000
    assert isinstance(query, str)

    #Choose the correct database
    if database == "reviewed":
        db=os.path.join(
                settings.BASE_DIR, 
                'pipeline', 'data', 'blast', 'clustercad_subunits_reviewed',
            )
    else:
        db=os.path.join(
                settings.BASE_DIR, 
                'pipeline', 'data', 'blast', 'clustercad_subunits_all',
            )

    # convert query to fasta format 
    queryStringIO = StringIO()
    queryRecord = SeqRecord(Seq(query, IUPAC.protein), id='query')
    SeqIO.write(queryRecord, queryStringIO, "fasta")
    queryFasta = queryStringIO.getvalue()
    queryStringIO.close()

    # run blast
    start = time()
    blastp_cline = NcbiblastpCommandline(
                                         db=db,
                                         evalue=evalue,
                                         outfmt=5,
                                         num_threads=2
                                         )
    result, stderr = blastp_cline(stdin=queryFasta)

    # parse blast output and delete files
    resultIO = StringIO(result)
    blast_record = NCBIXML.read(resultIO)
    resultIO.close()
    end = time()

    # iterate over record and generate output structure 
    alignments = []
    for alignment in blast_record.alignments:
        alignmentHit = alignment.title.split()[0].split('_')
        mibigAccession = alignmentHit[0]
        subunitId = alignmentHit[1]
        subunit = Subunit.objects.get(id=subunitId, cluster__mibigAccession=mibigAccession)
        hsps = []
        for hsp in alignment.hsps:
            domains = Domain.objects.filter(module__subunit=subunit, 
                                            stop__gte=hsp.sbjct_start,
                                            start__lte=hsp.sbjct_end).select_subclasses().order_by('start') 
            modules = list(set([domain.module for domain in domains]))
            modules = sorted(modules, key=lambda module: module.order)
            modules = [{'module': module, 'domains': list(domains.filter(module=module))} for module in modules]
            hsps.append({'hsp': hsp, 'modules': modules})
        alignments.append({'alignment': alignment, 'subunit': subunit, 'hsps': hsps})

    # if sortOutput=True, break apart HSPs and resort by bit order
    if sortOutput:
        individualHSPs = []
        for alignment in alignments:
            for hsp in alignment['hsps']:
                alignmentCopy = deepcopy(alignment)
                alignmentCopy['hsps'] = [hsp]
                individualHSPs.append(alignmentCopy)
        alignments = sorted(individualHSPs, key=lambda alignment: alignment['hsps'][0]['hsp'].bits, reverse=True)[0:max_target_seqs]

    return alignments, str(int(end-start))
