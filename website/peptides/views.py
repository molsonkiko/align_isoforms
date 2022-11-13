import json
import os
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.shortcuts import get_object_or_404, render
from django.urls import reverse
from django.views import generic

from .align_isoforms import get_isoforms, get_all_seqs, request_multi_alignment, get_protein as get_protein_json

from .models import Protein, Peptide, Isoform, Alignment

class IndexView(generic.ListView):
    template_name = "peptides/index.html"
    context_object_name = 'proteins'

    def get_queryset(self):
        '''get only the primary isoform of each protein'''
        return Protein.objects.exclude(acc_num__contains = '-')


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
            return Http404(
                reverse('peptides:proteins', args=(existing_prot.acc_num,))
            )
        except:
            pass
        prot = get_protein_json(acc_num)
        if not prot:
            return HttpResponse(
                "Could not find UniProt data for the accession number %s" % acc_num
            )
        prots = get_isoforms(prot)
        prot_seqs = get_all_seqs(prots)
        og_acc_num = acc_num
        og_prot = Protein(
            acc_num = og_acc_num, 
            sequence = prot_seqs[og_acc_num]
        )
        og_prot.save(force_insert=True)
        prot_list = og_acc_num
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
            Isoform.objects.create(
                prot_1 = Protein.objects.get(acc_num = og_acc_num),
                prot_2 = Protein.objects.get(acc_num = acc_num)
            )
        try:
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
    prot_objs = Protein.objects.filter(acc_num__in = acc_num_list)
    peptides =  (Peptide.objects
        .filter(prot__in = acc_num_list)
        .order_by('prot', 'location')
    )
    return render(
        request,
        'peptides/alignment.html',
        context = {
            'alignment': alignment,
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


def get_protein_json(request, acc_num: str):
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

def get_protein_json_schema(request):
    file_dir = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(file_dir, 'static/peptides/protein_json_schema.json')) as f:
        schema = f.read()
    return HttpResponse(
        schema,
        headers = {'content-type': 'application/json'}
    )