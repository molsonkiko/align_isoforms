import json
from pathlib import Path
import random
from django.test import TestCase

from .models import Protein, Peptide, Alignment, Isoform
from .sequence_chunkers import sequence_chunks, process_clustal_num
from .views import get_all_data_related_to_prot

CODE_DIR = Path(__file__).parent


class WebsiteTests(TestCase):
    maxDiff = 1000
    successfully_downloaded_P56856 = False
    def setUpTestData():
        '''Create a simple database with three isoforms,
        their peptides, and their multiple sequence
        alignment
        '''
        bluten_1 = Protein(
            sequence = 'MAWGKPRLFVCGTIK',
            acc_num = 'BLUTEN'
        )
        bluten_2 = Protein(
            sequence = 'MVTGKPRLTIK',
            acc_num = 'BLUTEN-2'
        )
        bluten_3 = Protein(
            sequence = 'CGTIR',
            acc_num = 'BLUTEN-3'
        )
        prots = [bluten_1, bluten_2, bluten_3]
        random.shuffle(prots)
        first_prot = prots[0]
        for prot in prots:
            prot.save()
        for other_prot in prots[1:]:
            Isoform.objects.create(
                prot_1 = first_prot,
                prot_2 = other_prot
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
        alignment = '''
CLUSTAL O(1.2.4) multiple sequence alignment

'''
        other_lines = [
            'BLUTEN-3      ----------CGTIR',
            'BLUTEN        MAWGKPRLFVCGTIK',
            'BLUTEN-2      MVTGKPRL----TIK'
        ]
        random.shuffle(other_lines)
        alignment += '\n'.join(other_lines)
        # the out-of-order proteins is realistic
        alignment += '\n' + ' '*14 + '*.:*****  ****.'
        Alignment.objects.create(
            prots = 'BLUTEN-3,BLUTEN,BLUTEN-2',
            alignment = alignment
        )
        # also download one real protein from UniProt and add its peptide
        Peptide.objects.create(
            prot = 'P56856',
            peptide = 'TSVFQYEGLWR'
        )
        __class__.successfully_downloaded_P56856 = get_all_data_related_to_prot('P56856')
        # finally create one-real-seeming protein with no peptides
        Protein.objects.create(acc_num = 'P56854', sequence = 'MMM')

    ### MODEL TESTS ###

    def test_isoforms_symmetric(self):
        '''verify that all three BLUTEN isoforms have the other two as isoforms
        '''
        bluten_1 = Protein.objects.get(acc_num = 'BLUTEN')
        bluten_2 = Protein.objects.get(acc_num = 'BLUTEN-2')
        bluten_3 = Protein.objects.get(acc_num = 'BLUTEN-3')
        self.assertEqual(set(x.acc_num for x in bluten_1.get_isoforms()), {'BLUTEN-2', 'BLUTEN-3'})
        self.assertEqual(set(x.acc_num for x in bluten_2.get_isoforms()), {'BLUTEN', 'BLUTEN-3'})
        self.assertEqual(set(x.acc_num for x in bluten_3.get_isoforms()), {'BLUTEN', 'BLUTEN-2'})

    def test_chunk_sequence(self):
        '''verify that sequence_chunkers.sequence_chunks
        works correctly'''
        bluten_1 = Protein.objects.get(acc_num = 'BLUTEN')
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
        self.assertTrue(self.successfully_downloaded_P56856)
        p56856 = Protein.objects.get(acc_num = 'P56856')
        p56856_isoforms = [x.acc_num for x in p56856.get_isoforms()]
        self.assertEqual(p56856_isoforms, ['P56856-2'])
        align = Alignment.objects.get(prots = 'P56856,P56856-2')
        self.assertIsNotNone(align)

    def test_chunk_alignment(self):
        bluten_3 = Protein.objects.get(acc_num = 'BLUTEN-3')
        peps = Peptide.objects.filter(prot__contains = 'BLUTEN')
        align = bluten_3.get_alignments()[0]
        chunks = process_clustal_num(align.alignment, peps, 9)
        second_chunk = chunks['chunks'][1]
        correct_chunk = [
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
        self.assertEqual(
            second_chunk,
            correct_chunk
        )

    def test_get_isoforms_in_acc_num_order(self):
        acc_nums = list(range(1, 15))
        random.shuffle(acc_nums)
        for ii in acc_nums:
            Protein.objects.create(
                acc_num = 'HUNDAR' if ii == 1 else 'HUNDAR-%i' % ii,
                sequence = 'MMM'
            )
        hundar = Protein.objects.get(acc_num = 'HUNDAR')
        acc_nums.remove(1)
        for ii in acc_nums:
            Isoform.objects.create(
                prot_1 = hundar,
                prot_2 = Protein.objects.get(acc_num = 'HUNDAR-%s' % ii),
            )
        hundar_isos = hundar.get_isoforms()
        self.assertEqual([x.isoform_num for x in hundar_isos],
            list(range(2, 15)))

    ### UI Tests ###

    def test_protein_page(self):
        response = self.client.get('/proteins/BLUTEN/?width=5')
        html = response.content.decode()
        isoform_html = '''
        <ol>
            <li> <a href = "/proteins/BLUTEN-2">Isoform BLUTEN-2</a> </li>
            <li> <a href = "/proteins/BLUTEN-3">Isoform BLUTEN-3</a> </li>
        </ol>
        '''
        self.assertInHTML(isoform_html, html)
        pep_spans = [
            '<span class="peptide 0" id = "0_0">MAW</span>',
            '<span class="peptide 6" id = "6_0">RLFV</span>',
            '<span class="peptide 6" id = "6_1">CG</span>',
            '<span class="peptide 11" id = "11_0">TI</span>',
        ]
        num_spans = [
            '<span class="right" style="left:55px">5</span>',
            '<span class="right" style="left:55px">10</span>',
            '<span class="right" style="left:55px">15</span>',
        ]
        for span in pep_spans + num_spans:
            with self.subTest(span = span, html = html):
                self.assertInHTML(span, html)
        non_pep_spans = ['GK', 'P', 'K']
        for span in non_pep_spans:
            with self.subTest(span = span, html = html):
                self.assertIn(f'-->{span}<!--', html)

    def test_protein_page_no_peps(self):
        response = self.client.get('/proteins/BLUTEN-3/')
        html = response.content.decode()
        self.assertIn('-->CGTIR<!--', html)

    def test_protein_json(self):
        response = self.client.get('/proteins/BLUTEN-2/json')
        blujson = json.loads(response.content.decode())
        alignment = Protein.objects.get(acc_num = 'BLUTEN').get_alignments()[0]
        correct_blujson = {
            'UniProt ID': 'BLUTEN-2',
            'sequence': 'MVTGKPRLTIK',
            'Isoform UniProt IDs': [
                'BLUTEN',
                'BLUTEN-3'
            ],
            'isoform alignments': [alignment.alignment],
            'mass spec peptides': [
                {
                    'sequence': 'PRLT',
                    'location': 5,
                }
            ]
        }
        self.assertEqual(blujson, correct_blujson)

    def test_get_protein(self):
        Peptide.objects.create(
            prot = 'Q9Y6Q5',
            peptide = 'VLFELTGR'
        )
        response = self.client.post('/get_protein/',
            {'acc_num': 'Q9Y6Q5'},
            follow=True # follow HTTP redirects; this is necessary
        )
        html = response.content.decode()
        # primary header
        self.assertInHTML(
            '<a href="/alignments/Q9Y6Q5,Q9Y6Q5-2">Alignment Q9Y6Q5,Q9Y6Q5-2</a>', 
            html
        )
        # sequence length header
        self.assertInHTML('<p>Sequence (length 423):</p>', html)
        # peptide link
        self.assertInHTML('<span class="sequence">VLFELTGR</span>', html)

    def test_json_schema(self):
        with open(CODE_DIR / 'static/peptides/protein_json_schema.json') as f:
            schema = json.load(f)
        response = self.client.get('/proteins/json_schema')
        got_schema = json.loads(response.content.decode())
        self.assertEqual(got_schema, schema)

    def test_alignments_page(self):
        response = self.client.get('/alignments/BLUTEN-3,BLUTEN,BLUTEN-2/')
        html = response.content.decode()
        pep_spans = [
            '<span class="peptide 0_BLUTEN" id="0_BLUTEN_0">MAW</span>',
            '<span class="peptide 6_BLUTEN" id="6_BLUTEN_0">RLFVCG</span>',
            '<span class="peptide 11_BLUTEN" id="11_BLUTEN_0">TI</span>',
            '<span class="peptide 5_BLUTEN-2" id="5_BLUTEN-2_0">PRL----T</span>'
        ]
        for span in pep_spans:
            with self.subTest(span = span, html = html):
                self.assertInHTML(span, html)
        self.assertInHTML('<span class="left-buffer">BLUTEN-3</span>', html)
        self.assertIn('-->*.:*****  ****.<!--', html)

    def test_download_alignment(self):
        disposition =  'attachment; filename = "isoforms of BLUTEN alignment.clustal_num"'
        content_type = 'text; charset = "utf-8"'
        response = self.client.get(
            '/download_alignment/BLUTEN-3,BLUTEN,BLUTEN-2/',
            headers = {
                'content-disposition': disposition,
                'content-type': content_type
            }
        )
        bluten_align = (Protein.objects
            .get(acc_num = 'BLUTEN')
            .get_alignments()
            [0].alignment
        )
        self.assertEqual(response.content.decode(), bluten_align)

    def test_peptide_csv_acc_num_like(self):
        response = self.client.get('/peptides/?acc_num_like=BLUTEN',
            headers = {
                'content-disposition': 'attachment; filename = "peptides of uniprot ids like BLUTEN.csv"',
                'content-type': 'text; charset = "utf-8"',
            }
        )
        lines = response.content.decode().split()
        correct_lines = [
            'uniprot_id,location,sequence',
            'BLUTEN,0,MAW',
            'BLUTEN,6,RLFVCG',
            'BLUTEN,11,GTI',
            'BLUTEN-2,5,PRLT'
        ]
        self.assertEqual(lines, correct_lines)

    def test_peptide_csv_acc_num_exact(self):
        response = self.client.get('/peptides/?acc_num=P56856',
            headers = {
                'content-disposition': 'attachment; filename = "peptides of uniprot id P56856.csv"',
                'content-type': 'text; charset = "utf-8"',
            }
        )
        lines = response.content.decode().split()
        correct_lines = [
            'uniprot_id,location,sequence',
            'P56856,40,TSVFQYEGLWR'
        ]
        self.assertEqual(lines, correct_lines)

    def test_peptide_csv_acc_num_exact_invalid_acc_num(self):
        response = self.client.get('/peptides/?acc_num=ZZZZZZZZZZZZ')
        html = response.content.decode()
        self.assertInHTML('Invalid accession number', html)

    def test_peptide_csv_acc_num_exact_not_in_db(self):
        response = self.client.get('/peptides/?acc_num=Q09996')
        html = response.content.decode()
        self.assertIn('No proteins match that query.', html)

    def test_peptide_csv_acc_num_exact_in_db_but_no_peps(self):
        response = self.client.get('/peptides/?acc_num=P56854')
        html = response.content.decode()
        self.assertIn(('A protein in the database has that accession number, '
                    'but no peptides in the database are associated with it.'), html)

    def test_peptide_csv_acc_num_like_in_db_but_no_peps(self):
        response = self.client.get('/peptides/?acc_num_like=P56854')
        html = response.content.decode()
        self.assertIn(('At least one protein in the database has that accession number, '
                    'but no peptides in the database are associated with them.'), html)

    def test_peptide_csv_acc_num_like_invalid_acc_num(self):
        response = self.client.get('/peptides/?acc_num_like=ZZZZZZZZZZZZ')
        html = response.content.decode()
        self.assertIn('Invalid accession number pattern', html)

    def test_peptide_csv_all(self):
        disposition = 'attachment; filename = "peptides.csv"'
        content_type = 'text; charset = "utf-8"'
        response = self.client.get('/peptides/',
            headers = {
                'content-type': content_type,
                'content-disposition': disposition
            }
        )
        lines = response.content.decode().split()
        correct_lines = [
            'uniprot_id,location,sequence',
            'BLUTEN,0,MAW',
            'BLUTEN,6,RLFVCG',
            'BLUTEN,11,GTI',
            'BLUTEN-2,5,PRLT',
            'P56856,40,TSVFQYEGLWR'
        ]
        self.assertEqual(lines, correct_lines)

    def test_index_default(self):
        response = self.client.get('/')
        html = response.content.decode()
        prot_list = '''
<ul>
    <li> <a href="proteins/BLUTEN">BLUTEN</a>: (length 15, 3 isoforms, 3 peptides) </li>
    <li> <a href="proteins/P56854">P56854</a>: (length 3, 1 isoform, 0 peptides) </li>
    <li> <a href="proteins/P56856">P56856</a>: (length 261, 2 isoforms, 1 peptide) </li>
</ul>'''
        self.assertInHTML(prot_list, html)

    def test_index_sort_by_len(self):
        response = self.client.get('/?orderby=len')
        html = response.content.decode()
        prot_list = '''
<ul>
    <li> <a href="proteins/P56854">P56854</a>: (length 3, 1 isoform, 0 peptides) </li>
    <li> <a href="proteins/BLUTEN">BLUTEN</a>: (length 15, 3 isoforms, 3 peptides) </li>
    <li> <a href="proteins/P56856">P56856</a>: (length 261, 2 isoforms, 1 peptide) </li>
</ul>'''
        self.assertInHTML(prot_list, html)

    def test_index_sort_by_num_isoforms(self):
        response = self.client.get('/?orderby=iso')
        html = response.content.decode()
        prot_list = '''
<ul>
    <li> <a href="proteins/P56854">P56854</a>: (length 3, 1 isoform, 0 peptides) </li>
    <li> <a href="proteins/P56856">P56856</a>: (length 261, 2 isoforms, 1 peptide) </li>
    <li> <a href="proteins/BLUTEN">BLUTEN</a>: (length 15, 3 isoforms, 3 peptides) </li>
</ul>'''
        self.assertInHTML(prot_list, html)
