# lib libraries
import time
# 3rd-party libraries
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests

# see https://rest.uniprot.org/docs/#/
BASE_QUERY = "https://rest.uniprot.org/uniprotkb/search?query=accession%3D"
def get_protein(acc_num: str) -> dict:
    '''Get the information in UniProt associated with accession number acc_num.
    API documentation: https://rest.uniprot.org/docs/#/uniprotkb/searchCursor
    '''
    resp = requests.get(BASE_QUERY + acc_num)
    resp.raise_for_status()
    return resp.json()['results']

def get_sequence(prot: dict) -> str:
    '''prot: JSON from the UniProt API for a protein
    Returns: the protein's sequence as a string
    '''
    return prot[0]['sequence']['value']

def get_isoform_ids(prot: dict) -> list:
    '''prot: JSON from the UniProt API for a protein
    Returns: a list of the UniProt accession nums for all isoforms
    '''
    comments = prot[0]['comments']
    out = set()
    for comment in comments:
        isoforms = comment.get('isoforms')
        if not isoforms:
            continue
        for isoform in isoforms:
            iso_ids = isoform.get('isoformIds')
            if not iso_ids:
                continue
            for iso_id in iso_ids:
                out.add(iso_id)
    return sorted(out)

def get_isoforms(prot: dict) -> dict:
    '''prot: JSON from the UniProt API for a protein
    Returns: A mapping of UniProt accession nums to the UniProt API JSON
    for all isoforms of the protein
    '''
    iso_ids = get_isoform_ids(prot)
    seqs = {}
    for id_ in iso_ids:
        try:
            seqs[id_] = get_protein(id_)
        except:
            continue
    return seqs

def get_all_prots(acc_num: str) -> dict:
    '''acc_num: The UniProt accession number of a protein
    Returns: a dict mapping the UniProt accession numbers for all distinct
    isoforms of a protein to the UniProt API JSON for that isoform.
    '''
    prot = get_protein(acc_num)
    prots = {acc_num: prot}
    seq = get_sequence(prot)
    isos = get_isoforms(prot)
    # remove the isoforms with the same sequence as the base acc num
    for iso_acc_num, iso in list(isos.items()):
        iso_seq = get_sequence(iso)
        if iso_acc_num != acc_num and iso_seq == seq:
            del isos[iso_acc_num]
    return {**prots, **isos}

def get_all_seqs(prots: dict) -> dict:
    '''prots: a dict mapping accession numbers to UniProt API JSON
    (returned by get_all_prots).
    Returns: a dict mapping those accession numbers to the sequences
    '''
    return {k: get_sequence(v) for k, v in prots.items()}

def to_fasta(seqs: dict) -> str:
    '''seqs: an {accession number -> UniProt API JSON} dict returned by
    get_all_prots.
    Returns: The accession numbers and sequences of those proteins in FASTA
    format.
    '''
    fasta = ''
    for acc_num, seq in seqs.items():
        seqrec = SeqRecord(Seq(seq))
        seqrec.id = acc_num
        fasta += seqrec.format('fasta')
    return fasta

# https://rest.uniprot.org/beta/docs/
WEBSITE_API = "https://rest.uniprot.org/beta"
def request_multi_alignment(seqs: dict) -> str:
    '''Send a request to the European Bioinformatics Institute
    for multiple alignment of several sequences
    '''
    fasta = to_fasta(seqs)
    r = requests.post(
        "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run", 
        data={
            "email": "lestimpe@gmail.com",
            "iterations": 1,
            "outfmt": "clustal_num",
            "order": "aligned",
            "sequence": fasta
           }
    )
    r.raise_for_status()
    job_id = r.text
    job_status = 'RUNNING'
    # ping the server every few seconds to see if the job is done
    while job_status == 'RUNNING':
        time.sleep(4)
        job_status_req = requests.get(
            f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}")
        job_status_req.raise_for_status()
        job_status = job_status_req.text
    # now that the job is done, get the alignment
    resp = requests.get(
        f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-clustal_num")
    resp.raise_for_status()
    return resp.text

def align_isoforms(acc_num: str, **kwargs) -> str:
    '''get all isoforms of the protein with accession number acc_num,
    and return a sequence alignment.
    Can also pass in keyword arguments to set parameters for the alignment,
    if there are only two isoforms to align'''
    prots = get_all_prots(acc_num)
    seqs = get_all_seqs(prots)
    if len(seqs) < 2:
        print(f"The protein with UniProt accession number {acc_num} has only one isoform")
        (acc1, seq1) = list(seqs.items())[0]
        return seq1
    return request_multi_alignment(seqs)