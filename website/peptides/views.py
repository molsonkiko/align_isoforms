from django.http import HttpResponseRedirect
from django.shortcuts import get_object_or_404, render
from django.urls import reverse
from django.views import generic

from .align_isoforms import get_all_prots, get_all_seqs, request_multi_alignment

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
            return HttpResponseRedirect(
                reverse('peptides:proteins', args=(existing_prot.acc_num,))
            )
        except:
            pass
        prots = get_all_prots(acc_num)
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
    peptides =  Peptide.objects.filter(prot__in = acc_num_list)
    return render(
        request,
        'peptides/alignment.html',
        context = {
            'alignment': alignment,
            'proteins': prot_objs,
            'peptides': peptides,
        }
    )