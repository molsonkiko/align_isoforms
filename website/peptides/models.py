import re
from django.contrib import admin
from django.db import models

class BaseModel(models.Model):
    class Meta:
        app_label = 'peptides'
        abstract = True


is_acc_num = re.compile(
    ("(?:[OPQ][0-9][A-Z0-9]{3}[0-9]" # first style, e.g., P56856
    "|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})" # second style, e.g., A2BC19, A0A023GPI8
    "(?:-([1-9][0-9]*))?") # optional isoform number, e.g., -2, -10
).fullmatch

def isoform_num(acc_num: str) -> int:
    acc_num_split = acc_num.split('-')
    if len(acc_num_split) == 1:
        return 1
    return int(acc_num_split[1])


class Protein(BaseModel):
    prot_id = models.AutoField(primary_key=True)
    acc_num = models.CharField(max_length=15, unique=True)
    sequence = models.CharField(max_length=45_000)
    isoform_num = models.IntegerField(default = 1)
    # longest protein ever found is titin, with ~35,000 aa

    class Meta:
        indexes = [
            models.Index(fields = ['acc_num'], name = 'acc_num_idx')
        ]

    def __str__(self) -> str:
        return self.acc_num

    def __repr__(self) -> str:
        return 'Protein(%s)' % self.acc_num

    def save(self, *args, **kwargs):
        '''when protein saved, update the location of each associated 
        peptide to reflect its location in the protein.
        Also save the isoform number of the protein 
        (e.g., BLUTEN-1 has isoform # 1, BLUTEN-3 has isoform # 3)
        '''
        self.isoform_num = isoform_num(self.acc_num)
        super().save(*args, **kwargs)
        for peptide in Peptide.objects.filter(prot = self.acc_num):
            peptide.save(force_update=True)

    def get_isoforms(self):
        '''Return all protein objects that are isoforms of self
        '''
        # if this is the primary isoform, this will get all isoforms
        isos_2 = Isoform.objects.filter(prot_1 = self)
        if not isos_2:
            # this is a secondary isoform, or has no isoforms
            isos_1 = Isoform.objects.filter(prot_2 = self)
            if not isos_1:
                return [] # no isoforms
            # get the primary isoform, then get all isoforms of the primary
            primary_iso = isos_1.values('prot_1')[0]['prot_1']
            other_isos = Isoform.objects.filter(prot_1 = primary_iso)
            all_isos = isos_2.union(other_isos)
        else:
            all_isos = isos_2
        acc_nums = set()
        for iso in all_isos:
            acc_nums.add(iso.prot_1.acc_num)
            acc_nums.add(iso.prot_2.acc_num)
        acc_nums.remove(self.acc_num)
        return (Protein.objects.filter(acc_num__in = acc_nums)
            .order_by('isoform_num'))


    def get_alignments(self):
        '''return all alignment objects where self is one of the aligned
        proteins
        '''
        # equivalent to SELECT * FROM Alignment WHERE prots LIKE '%{acc_num}%'
        return Alignment.objects.filter(prots__contains = self.acc_num)

    def get_peptides(self):
        '''Return all peptides (in the Peptide table) associated with self
        '''
        return (Peptide.objects
            .filter(prot = self.acc_num)
            .order_by('location')
        )


class Isoform(BaseModel):
    pkey = models.AutoField(primary_key=True)
    prot_1 = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='prot_1')
    prot_2 = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='prot_2')

    def __str__(self) -> str:
        return "Isoform(%s, %s)" % (self.prot_1, self.prot_2)

    __repr__ = __str__


class Alignment(BaseModel):
    prots = models.CharField(max_length=300, primary_key=True)
    alignment = models.CharField(max_length=720_000)


class Peptide(BaseModel):
    pkey = models.AutoField(primary_key=True)
    prot = models.CharField(max_length=15)
    peptide = models.CharField(max_length=10_000)
    location = models.IntegerField(default=-1)

    def __str__(self):
        if len(self.peptide) > 20:
            prot_abbrev = self.peptide[:17] + '...'
        else: prot_abbrev = self.peptide
        return 'Peptide(%s, %i, %s)' % (self.prot, self.location, prot_abbrev)

    __repr__ = __str__

    def save(self, *args, **kwargs):
        '''Find the protein with this peptide's accession number, and if that
        protein is in the database, set the location field to the location
        of this peptide in that protein's sequence.
        '''
        corresponding_protein = Protein.objects.filter(acc_num = self.prot)
        if len(corresponding_protein) == 0:
            super().save(*args, **kwargs)
            return
        try:
            self.location = corresponding_protein[0].sequence.index(self.peptide)
        except ValueError: # raised by .index if thing not in list
            pass
        super().save(*args, **kwargs)

    @admin.display
    def peptide_preview(self):
        if len(self.peptide) < 20:
            return self.peptide
        return self.peptide[:17] + '...'
    
    class Meta:
        indexes = [
            models.Index(fields = ['prot'], name = 'prot_idx')
        ]