from django.db import models

class BaseModel(models.Model):
    class Meta:
        app_label = 'peptides'
        abstract = True


class Protein(BaseModel):
    prot_id = models.AutoField(primary_key=True)
    acc_num = models.CharField(max_length=15, unique=True)
    sequence = models.CharField(max_length=45_000)
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
        peptide to reflect its location in the protein
        '''
        super().save(*args, **kwargs)
        peptides = Peptide.objects.filter(prot = self.acc_num)
        for peptide in peptides:
            pep = peptide.peptide
            peptide.location = self.sequence.index(pep)
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
            .order_by('acc_num'))


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
    prot_1 = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='prot_1')
    prot_2 = models.ForeignKey(Protein, on_delete=models.CASCADE, related_name='prot_2')

    def __str__(self) -> str:
        return "Isoform(%s, %s)" % (self.prot_1, self.prot_2)

    __repr__ = __str__


class Alignment(BaseModel):
    prots = models.CharField(max_length=250, primary_key=True)
    alignment = models.CharField(max_length=720_000)


class Peptide(BaseModel):
    prot = models.CharField(max_length=15)
    peptide = models.CharField(max_length=10_000)
    location = models.IntegerField(default=-1)

    def __str__(self):
        if len(self.peptide) > 20:
            prot_abbrev = self.peptide[:17] + '...'
        else: prot_abbrev = self.peptide
        return 'Peptide(%s, %i, %s)' % (self.prot, self.location, prot_abbrev)

    __repr__ = __str__
    
    class Meta:
        indexes = [
            models.Index(fields = ['prot'], name = 'prot_idx')
        ]