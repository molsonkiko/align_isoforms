import re
from django.contrib import admin
from django.http import HttpResponseRedirect
from django.shortcuts import render
from django.urls import path

from .models import Protein, Alignment, Isoform, Peptide

admin.site.register(Protein)
admin.site.register(Alignment)
admin.site.register(Isoform)


class PeptideAdmin(admin.ModelAdmin):
    list_display = ('prot', 'location', 'peptide_preview')

    def get_urls(self):
        urls = super().get_urls()
        my_urls = [
            path('from_csv/', admin.site.admin_view(peptides_from_csv))
        ]
        return my_urls + urls

is_peptide_str = re.compile('[ACDEFGHIKLMNPQRSTVWY]+').fullmatch

def peptides_from_csv(request, *args, **kwargs):
    csv_file = request.FILES.get('csv_file')
    # this is an UploadedFile object
    # see https://docs.djangoproject.com/en/4.1/ref/files/uploads/#django.core.files.uploadedfile.UploadedFile
    if not csv_file:
        return render(
            request, 
            'admin/peptides_from_csv.html',
            context={
                'error_message': 'No CSV file was selected'
            }
        )
    sep = request.POST.get('sep', ',')
    if sep.lower() == 'tab' or sep == '\\t':
        sep = '\t'
    csv_file = csv_file.read().decode()
    peps = []
    for line in re.split('\r?\n|\r', csv_file):
        fields = line.split(sep)
        if len(fields) < 2:
            return render(
                request,
                'admin/peptides_from_csv.html',
                context={
                    'error_message': 'CSV file must have two or three columns delimited by %r' % sep
                }
            )
        seq_field = fields[1]
        fields[1] = seq_field.upper()
        if not is_peptide_str(seq_field):
            # we assume that an all-caps series of amino acid letters is just
            # an amino acid sequence and not a header
            continue
        if len(fields) == 3:
            try:
                loc = int(fields[2])
            except:
                continue
            pep = Peptide(prot = fields[0], peptide = fields[1], location = loc)
        else:
            pep = Peptide(prot = fields[0], peptide = fields[1])
        peps.append(pep)
    Peptide.objects.bulk_create(peps)
    return HttpResponseRedirect('/admin/peptides/peptide')


admin.site.register(Peptide, PeptideAdmin)