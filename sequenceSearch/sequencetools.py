import os
from textwrap import wrap
from tempfile import mkstemp
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from pks.models import Subunit, Domain
from model_utils.managers import InheritanceManager

def blast(query, db="/clusterCAD/pipeline/data/blast/clustercad_subunits", evalue=10.0, max_target_seqs=10):
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

    # write query to a temporary fasta file
    queryFasta = string2Fasta("query", query)
    queryFile = mkstemp(text=True)
    f = os.fdopen(queryFile[0], "w")
    f.write(queryFasta)
    f.close()

    # run blastp
    outputfile = mkstemp(text=True)
    blastp_cline = NcbiblastpCommandline(
                                         query=queryFile[1],
                                         db=db,
                                         evalue=evalue,
                                         outfmt=5,
                                         out=outputfile[1],
                                         max_target_seqs=max_target_seqs)
    stdout, stderr = blastp_cline()

    # parse blastp output and delete files
    result_handle = open(outputfile[1])
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()
    os.remove(queryFile[1])
    os.remove(outputfile[1])

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
                                            start__lte=hsp.sbjct_end).select_subclasses() 
            modules = list(set([domain.module for domain in domains]))
            modules = sorted(modules, key=lambda module: module.order)
            modules = [{'module': module, 'domains': list(domains.filter(module=module))} for module in modules]
            hsps.append({'hsp': hsp, 'modules': modules})
        alignments.append({'alignment': alignment, 'subunit': subunit, 'hsps': hsps})

    return alignments

def string2Fasta(headers, strings):
    # convert a set of strings and headers to a fasta file

    if isinstance(headers, str):
        headers = [headers]
    if isinstance(strings, str):
        strings = [strings]
    assert len(headers) == len(strings)

    zipped = zip(['>' + x for x in headers], strings)
    lines = ['\n'.join(wrap(x)) + '\n' for sublist in zipped for x in sublist]
    return ''.join(lines)
