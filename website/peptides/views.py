import json
import os
import re
from django.http import HttpResponseRedirect, HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from django.urls import reverse
# from django.views.decorators.cache import never_cache

from .align_isoforms import get_isoforms, get_all_seqs, request_multi_alignment, get_protein as uniprot_json, align_isoforms
from .models import Protein, Peptide, Isoform, Alignment, is_acc_num
from .sequence_chunkers import sequence_chunks, process_clustal_num

CODE_DIR = os.path.join(os.path.dirname(__file__))

def index_view(request):
    '''Show only the primary isoforms of proteins in the database.
    Optionally allow to order by length, by number of isoforms,
    or alphabetically.
    '''
    proteins = Protein.objects.filter(isoform_num = 1)
    data = [
        {
            'acc_num': prot.acc_num,
            'npeps': len(prot.get_peptides()),
            'lenseq': len(prot.sequence),
            'n_isoforms': len(prot.get_isoforms()) + 1,
        }
        for prot in proteins
    ]
    orderby = request.GET.get('orderby', 'alpha')
    sort_reverse = orderby[0] == '-'
    if sort_reverse:
        orderby = orderby[1:]
    orderby = {
        'alpha': 'acc_num',
        'iso': 'n_isoforms',
        'len': 'lenseq',
        'npeps': 'npeps',
    }[orderby] # sort by accession number unless user says otherwise
    return render(
        request,
        'peptides/index.html',
        context = {
            'prot_data': data, 
            'orderby': orderby, 
            'sort_reverse': "true" if sort_reverse else "false"
        }
        # sorting is done server-side; see script at top of index.html
    )


def about_view(request):
    return render(request, 'peptides/about.html')


def protein_view(request, acc_num: str):
    width = int(request.GET.get('width', 120))
    num_offset = 10 + 9 * width
    if acc_num.endswith('-1'):
        acc_num = acc_num[:-2]
        # the primary isoform may end with '-1', but it doesn't
        # matter and it's really just confusing
    prot = get_object_or_404(Protein, acc_num = acc_num)
    alignments = prot.get_alignments()
    isoforms = prot.get_isoforms()
    peptides = prot.get_peptides()
    chunks = sequence_chunks(prot.sequence, peptides, width)
    annotated_chunks = []
    chunk_end = 0
    for chunk in chunks:
        chunk_end += sum(len(p['seq']) for p in chunk)
        annotated_chunks.append(
            {'chunk': chunk, 'chunk_end': chunk_end}
        )
    return render(
        request,
        template_name='peptides/protein.html',
        context = {
            'protein': prot,
            'alignments': alignments,
            'isoforms': isoforms,
            'peptides': peptides,
            'sequence_chunks': annotated_chunks,
            'seq_len': len(prot.sequence),
            'num_offset': num_offset,
        }
    )


def get_all_data_related_to_prot(acc_num: str) -> bool:
    '''Get all isoforms of protein with UniProt accession number acc_num,
    get a multiple sequence alignment of those isoforms,
    and for each peptide in the database that belongs to one of the isoforms,
    get its location in the sequence of that isoform.
    '''
    # get uniprot data for the protein and save it to the database
    prot = uniprot_json(acc_num)
    if not prot:
        return False
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
    return True


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
        except: # no existing prot, so get the data
            pass
        try:
            data = get_all_data_related_to_prot(acc_num)
            if not data:
                return HttpResponse(
                    "Could not find UniProt data for the accession number " + acc_num
                )
        except:
            response = HttpResponse(
                ("The server had an error while retrieving data on protein %s "
                "(or its isoforms) from UniProt.") % acc_num
            )
            response.status_code = 500
            return response
        return HttpResponseRedirect(
            reverse('peptides:proteins', args=(acc_num,))
        )


def alignments_view(request, acc_nums: str):
    width = int(request.GET.get('width', 60))
    num_offset = 120 + 9 * width
    alignment = get_object_or_404(Alignment, pk=acc_nums)
    acc_num_list = acc_nums.split(',')
    prot_objs = (Protein.objects
        .filter(acc_num__in = acc_num_list)
        .order_by('isoform_num')
    )
    peptides =  (Peptide.objects
        .filter(prot__in = acc_num_list)
        .order_by('prot', 'location')
    )
    alignment_pieces = process_clustal_num(alignment.alignment, peptides, width)
    return render(
        request,
        'peptides/alignment.html',
        context = {
            'prots': alignment.prots,
            'alignment_pieces': alignment_pieces,
            'proteins': prot_objs,
            'peptides': peptides,
            'num_offset': num_offset,
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
    return JsonResponse(
        {
            'UniProt ID': prot.acc_num,
            'sequence': prot.sequence,
            'Isoform UniProt IDs': isoform_ids,
            'isoform alignments': alignments,
            'mass spec peptides': peptides, 
        }
    )


def protein_json_schema(request):
    with open(os.path.join(CODE_DIR, 'static/peptides/protein_json_schema.json')) as f:
        schema = json.load(f)
    return JsonResponse(schema)


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
        if not re.fullmatch(r'[A-Z\d-]{6,10}', acc_num_like):
            return HttpResponse('Invalid accession number pattern')
        filename =  'peptides of uniprot ids like %s.csv' % acc_num_like
        peptides = (Peptide.objects
            .filter(prot__contains = acc_num_like)
            .order_by('prot', 'location')
        )
    elif acc_num:
        if not is_acc_num(acc_num):
            return HttpResponse('Invalid accession number')
        filename = 'peptides of uniprot id %s.csv' % acc_num
        peptides = (Peptide.objects
            .filter(prot = acc_num)
            .order_by('prot', 'location')
        )
    else:
        filename = 'peptides.csv'
        peptides = Peptide.objects.order_by('prot', 'location')
        if len(peptides) == 0:
            return HttpResponse('There are no peptides currently in the database.')
    if len(peptides) == 0:
        prots_with_acc_num = None
        if acc_num:
            prots_with_acc_num = Protein.objects.filter(acc_num = acc_num)
            if prots_with_acc_num:
                return HttpResponse(
                    ('A protein in the database has that accession number, '
                    'but no peptides in the database are associated with it.')
                )
        elif acc_num_like:
            prots_with_acc_num = Protein.objects.filter(acc_num__contains = acc_num)
            if prots_with_acc_num:
                return HttpResponse(
                    ('At least one protein in the database has that accession number, '
                    'but no peptides in the database are associated with them.')
                )
        return HttpResponse('No proteins match that query.')
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


def request_alignment(request):
    '''Sometimes the EBI computer won't return an alignment when the user gets the
    data for a protein and its isoforms.
    This allows the user to resubmit a request for data from the EBI.
    When the request succeeds, or if the protein's alignment was already in 
    the database, redirect to the alignment page.
    '''
    acc_num = request.POST.get('acc_num')
    if not acc_num:
        return HttpResponse("Must supply an accession number when requesting an alignment.")
    if not is_acc_num(acc_num):
        return HttpResponse(f"Can't request an alignment for '{acc_num}' because it's not a valid accession number.")
    prot = Protein.objects.get(acc_num = acc_num)
    existing_alignment = prot.get_alignments()
    isoforms = prot.get_isoforms()
    prot_list = ','.join([acc_num] + [x.acc_num for x in isoforms])
    if existing_alignment:
        return HttpResponseRedirect('/alignments/' + prot_list)
    try:
        alignment = align_isoforms(acc_num)
    except Exception as ex:
        return HttpResponse("While requesting alignment, got the following error: " + str(ex))
    if not alignment:
        return HttpResponse("The EBI did not return an alignment for protein '" + acc_num + "'. Try again later.")
    Alignment.objects.create(prots=prot_list, alignment=alignment)
    return HttpResponseRedirect("/alignments/" + prot_list)