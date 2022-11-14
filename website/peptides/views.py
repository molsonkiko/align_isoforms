import json
import os
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.shortcuts import get_object_or_404, render
from django.urls import reverse
from django.views import generic

from .align_isoforms import get_isoforms, get_all_seqs, request_multi_alignment, get_protein as uniprot_json
from .models import Protein, Peptide, Isoform, Alignment
from .sequence_chunkers import sequence_chunks, process_clustal_num

def index_view(request):
    '''Show only the primary isoforms of proteins in the database.
    Optionally allow to order by length, by number of isoforms,
    or alphabetically.
    '''
    proteins = Protein.objects.exclude(acc_num__contains = '-')
    orderby = request.GET.get('orderby')
    if orderby:
        if orderby == 'alpha':
            proteins = sorted(proteins, key = lambda x: x.acc_num)
        elif orderby[:3] == 'iso':
            proteins = sorted(proteins, key = lambda x: len(x.get_isoforms()))
        elif orderby == 'len':
            proteins = sorted(proteins, key = lambda x: len(x.sequence))
    return render(
        request,
        'peptides/index.html',
        context = {'proteins': proteins}
    )

def protein_view(request, acc_num: str):
    if acc_num.endswith('-1'):
        acc_num = acc_num[:-2]
        # the primary isoform may end with '-1', but it doesn't
        # matter and it's really just confusing
    prot = get_object_or_404(Protein, acc_num = acc_num)
    alignments = prot.get_alignments()
    isoforms = prot.get_isoforms()
    peptides = prot.get_peptides()
    return render(
        request,
        template_name='peptides/protein.html',
        context = {
            'protein': prot,
            'alignments': alignments,
            'isoforms': isoforms,
            'peptides': peptides,
            'sequence_chunks': sequence_chunks(prot.sequence, peptides),
            'seq_len': len(prot.sequence),
        }
    )


def get_protein(request):
    try:
        acc_num = request.POST['acc_num']
    except:
        return render(request, 'peptides/get_protein.html')
    else:
        if acc_num.endswith('-1'):
            acc_num = acc_num[:-2]
        try:
            existing_prot = Protein.objects.get(acc_num = acc_num)
            return HttpResponseRedirect(
                reverse('peptides:proteins', args=(existing_prot.acc_num,))
            )
        except:
            pass
        # get uniprot data for the protein and save it to the database
        prot = uniprot_json(acc_num)
        if not prot:
            return HttpResponse(
                "Could not find UniProt data for the accession number %s" % acc_num
            )
        og_acc_num = acc_num
        # now get all the uniprot data for all the isoforms of the protein 
        prots = get_isoforms(prot)
        prot_seqs = get_all_seqs(prots)
        og_prot = Protein(
            acc_num = og_acc_num, 
            sequence = prot_seqs[og_acc_num]
        )
        og_prot.save(force_insert=True)
        prot_list = og_acc_num
        # add each isoform to the database
        for acc_num, seq in prot_seqs.items():
            if acc_num == og_acc_num:
                continue
            if acc_num.endswith('-1'):
                acc_num = acc_num[:-2]
            prot = Protein(
                acc_num = acc_num, 
                sequence = seq
            )
            prot_list += ',' + acc_num
            prot.save(force_insert=True)
            # also add a relationschip between og_prot and new prot
            Isoform.objects.create(
                prot_1 = Protein.objects.get(acc_num = og_acc_num),
                prot_2 = Protein.objects.get(acc_num = acc_num)
            )
        try:
            # finally, get a multiple sequence alignment of all the isoforms
            alignment = request_multi_alignment(prot_seqs)
            Alignment.objects.create(prots = prot_list, alignment = alignment)
        except:
            pass
        return HttpResponseRedirect(
            reverse('peptides:proteins', args=(og_prot.acc_num,))
        )


def alignments_view(request, acc_nums: str):
    alignment = get_object_or_404(Alignment, pk=acc_nums)
    acc_num_list = acc_nums.split(',')
    prot_objs = (Protein.objects
        .filter(acc_num__in = acc_num_list)
        .order_by('acc_num')
    )
    peptides =  (Peptide.objects
        .filter(prot__in = acc_num_list)
        .order_by('prot', 'location')
    )
    alignment_pieces = process_clustal_num(alignment.alignment, peptides)
    return render(
        request,
        'peptides/alignment.html',
        context = {
            'prots': alignment.prots,
            'alignment_pieces': alignment_pieces,
            'proteins': prot_objs,
            'peptides': peptides,
        }
    )


def download_alignment(request, prots: str):
    alignment = Alignment.objects.get(pk = prots)
    primary_acc_num = alignment.prots.split(',')[0]
    try:
        dash_index = primary_acc_num.index('-')
        # all isoforms have the same alignment as the primary isoform
        primary_acc_num = primary_acc_num[:dash_index]
    except: pass
    # this disposition indicates that it should be downloaded with this default file name
    disposition =  'attachment; filename = "isoforms of %s alignment.clustal_num"' % primary_acc_num
    # this content-type indicates that this is a utf-8 text file
    content_type = 'text; charset = "utf-8"'
    return HttpResponse(
        alignment.alignment,
        headers = {'content-disposition': disposition, 'content-type': content_type},
    )


def protein_json(request, acc_num: str):
    prot = Protein.objects.get(acc_num = acc_num)
    isoform_ids = [iso.acc_num for iso in prot.get_isoforms()]
    alignments = [al.alignment for al in prot.get_alignments()]
    peptides = [{'location': pep.location, 'sequence': pep.peptide} 
        for pep in prot.get_peptides()]
    return HttpResponse(
        json.dumps({
            'UniProt ID': prot.acc_num,
            'sequence': prot.sequence,
            'Isoform UniProt IDs': isoform_ids,
            'isoform alignments': alignments,
            'mass spec peptides': peptides, 
        }),
        headers = {
            'content-type': 'application/json',
        }
    )


def protein_json_schema(request):
    file_dir = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(file_dir, 'static/peptides/protein_json_schema.json')) as f:
        schema = f.read()
    return HttpResponse(
        schema,
        headers = {'content-type': 'application/json'}
    )


def peptides_csv(request):
    '''Return a csv containing responsive accession numbers.
    Query must be of the form 'acc_num_like=<partial accession number>'
        or 'acc_num=<exact accession number>', or no query at all.
    If no query, return a csv returning all peptides.
    If acc_num query, return only the peptides for that accession number.
    If acc_num_like query, return peptides where the accession number
        contains the acc_num_like substring.
    '''
    acc_num_like = request.GET.get('acc_num_like', '')
    acc_num = request.GET.get('acc_num', '')
    if acc_num_like:
        filename =  'peptides of uniprot ids like %s.csv' % acc_num_like
        peptides = (Peptide.objects
            .filter(prot__contains = acc_num_like)
            .order_by('prot', 'location')
        )
    elif acc_num:
        filename = 'peptides of uniprot id %s.csv' % acc_num
        peptides = (Peptide.objects
            .filter(prot = acc_num)
            .order_by('prot', 'location')
        )
    else:
        filename = 'peptides.csv'
        peptides = Peptide.objects.order_by('prot', 'location')
    csv = 'uniprot_id,location,sequence\n'
    for pep in peptides:
        csv += '%s,%s,%s\n' % (pep.prot, pep.location, pep.peptide)
    disposition = 'attachment; filename = "%s"' % filename
    content_type = 'text; charset = "utf-8"'
    return HttpResponse(
        csv,
        headers = {'content-disposition': disposition, 'content-type': content_type},
    )


def site_map(request):
    return render(request, 'peptides/site_map.html')