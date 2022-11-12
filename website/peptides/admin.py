from django.contrib import admin

from .models import Protein, Alignment, Isoform, Peptide

admin.site.register(Protein)
admin.site.register(Alignment)
admin.site.register(Isoform)
admin.site.register(Peptide)