import json
from django.test import TestCase

from .models import Protein, Peptide, Alignment, Isoform
from .views import get_all_data_related_to_prot

from .sequence_chunkers import sequence_chunks, process_clustal_num

class ProteinTests(TestCase):
    def setUpTestData():
        '''Create a simple database with three isoforms,
        their peptides, and their multiple sequence
        alignment
        '''
        bluten_1 = Protein(
            sequence = 'MAWGKPRLFVCGTIK',
            acc_num = 'BLUTEN'
        )
        bluten_1.save()
        bluten_2 = Protein(
            sequence = 'MVTGKPRLTIK',
            acc_num = 'BLUTEN-2'
        )
        bluten_2.save()
        bluten_3 = Protein(
            sequence = 'CGTIR',
            acc_num = 'BLUTEN-3'
        )
        bluten_3.save()
        Isoform.objects.create(
            prot_1 = bluten_1,
            prot_2 = bluten_2
        )
        Isoform.objects.create(
            prot_1 = bluten_1,
            prot_2 = bluten_3
        )
        Peptide.objects.create(
            prot = 'BLUTEN',
            peptide = 'MAW',
            location = 0
        )
        Peptide.objects.create(
            prot = 'BLUTEN',
            peptide = 'RLFVCG',
            location = 6
        )
        Peptide.objects.create(
            prot = 'BLUTEN',
            peptide = 'GTI',
            location = 11
        )
        Peptide.objects.create(
            prot = 'BLUTEN-2',
            peptide = 'PRLT',
            location = 5
        )
        Alignment.objects.create(
            prots = 'BLUTEN,BLUTEN-2,BLUTEN-3',
            alignment = '''
CLUSTAL O(1.2.4) multiple sequence alignment

BLUTEN        MAWGKPRLFVCGTIK
BLUTEN-2      MVTGKPRL----TIK
BLUTEN-3      ----------CGTIR
              *.:*****  ****.'''
        )

    def test_isoforms_symmetric(self):
        '''verify that all three BLUTEN isoforms have the other two as isoforms
        '''
        bluten_1 = Protein.objects.get(pk = 1)
        bluten_2 = Protein.objects.get(pk = 2)
        bluten_3 = Protein.objects.get(pk = 3)
        self.assertEqual(set(x.acc_num for x in bluten_1.get_isoforms()), {'BLUTEN-2', 'BLUTEN-3'})
        self.assertEqual(set(x.acc_num for x in bluten_2.get_isoforms()), {'BLUTEN', 'BLUTEN-3'})
        self.assertEqual(set(x.acc_num for x in bluten_3.get_isoforms()), {'BLUTEN', 'BLUTEN-2'})

    def test_chunk_sequence(self):
        '''verify that sequence_chunkers.sequence_chunks
        works correctly'''
        bluten_1 = Protein.objects.get(pk = 1)
        peps = bluten_1.get_peptides()
        chunks = sequence_chunks(bluten_1.sequence, peps, 9)
        self.assertEqual(
            chunks,
            [
                [
                    {
                        'seq': '',
                        'is_pep': False
                    },
                    {
                        'seq': 'MAW', 
                        'is_pep': True,
                        'loc': 0,
                        'pep_num': 0
                    },
                    {
                        'seq': 'GKP',
                        'is_pep': False
                    },
                    {
                        'seq': 'RLF',
                        'is_pep': True,
                        'loc': 6,
                        'pep_num': 0
                    }
                ],
                [
                    {
                        'seq': 'VCG',
                        'is_pep': True,
                        'loc': 6,
                        'pep_num': 1
                    },
                    {
                        'seq': 'TI',
                        'is_pep': True,
                        'loc': 11,
                        'pep_num': 0
                    },
                    {
                        'seq': 'K',
                        'is_pep': False
                    }
                ]
            ]
        )

    def test_deleting_protein_deletes_isoforms(self):
        bluten = Protein.objects.get(acc_num = 'BLUTEN')
        bluten.delete()
        isoforms = Isoform.objects.filter(
            prot_1 = bluten
        ).union(Isoform.objects.filter(
            prot_2 = bluten
        ))
        self.assertEqual(len(isoforms), 0)

    def test_add_then_delete_then_add_again(self):
        foobar = Protein(
            sequence = 'FOOBARKL',
            acc_num = 'FOOBAR'
        )
        foobar.save()
        foobar_2 = Protein(
            sequence = 'FOARKL',
            acc_num = 'FOOBAR-2'
        )
        foobar_2.save()
        Isoform.objects.create(
            prot_1 = foobar,
            prot_2 = foobar_2
        )
        foobar.delete()
        self.assertEqual(list(foobar_2.get_isoforms()), [])
        foobar_2.delete()
        foobar = Protein(
            sequence = 'FOOBARKL',
            acc_num = 'FOOBAR'
        )
        foobar.save()
        foobar_2 = Protein(
            sequence = 'FOARKL',
            acc_num = 'FOOBAR-2'
        )
        foobar_2.save()
        Isoform.objects.create(
            prot_1 = foobar,
            prot_2 = foobar_2
        )
    
    def test_query_uniprot(self):
        self.assertTrue(get_all_data_related_to_prot('P56856'))
        p56856 = Protein.objects.get(acc_num = 'P56856')
        p56856_isoforms = [x.acc_num for x in p56856.get_isoforms()]
        self.assertEqual(p56856_isoforms, ['P56856-2'])
        align = Alignment.objects.get(prots = 'P56856,P56856-2')
        self.assertIsNotNone(align)

    def test_chunk_alignment(self):
        bluten_3 = Protein.objects.get(pk = 3)
        peps = Peptide.objects.filter(prot__contains = 'BLUTEN')
        align = bluten_3.get_alignments()[0]
        chunks = process_clustal_num(align.alignment, peps, 9)
        second_chunk = chunks['chunks'][1]
        # print(json.dumps(second_chunk, indent=4))
        self.assertEqual(
            second_chunk,
            [
                {
                    'acc_num': 'BLUTEN',
                    'chunk': [
                        {
                            'seq': 'VCG',
                            'is_pep': True,
                            'loc': 6,
                            'pep_num': 1
                        },
                        {
                            'seq': 'TI',
                            'is_pep': True,
                            'loc': 11,
                            'pep_num': 0
                        },
                        {
                            'seq': 'K',
                            'is_pep': False
                        }
                    ],
                    'chunk_end': 15,
                    'is_prot_chunk': True
                },
                {
                    'acc_num': 'BLUTEN-2',
                    'chunk': [
                        {
                            'seq': '---T',
                            'is_pep': True,
                            'loc': 5,
                            'pep_num': 1
                        },
                        {
                            'seq': 'IK',
                            'is_pep': False
                        }
                    ],
                    'chunk_end': 11,
                    'is_prot_chunk': True
                },
                {
                    'acc_num': 'BLUTEN-3',
                    'chunk': [
                        {
                            'seq': '-CGTIR',
                            'is_pep': False
                        }
                    ],
                    'chunk_end': 5,
                    'is_prot_chunk': True
                },
                {
                    'acc_num': '',
                    'chunk': [
                        {
                            'seq': ' ****.',
                            'is_pep': False
                        }
                    ],
                    'chunk_end': 15,
                    'is_prot_chunk': False
                }
            ]
        )