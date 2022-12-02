# Change Log
All [notable changes](#0112---2022-12-01) to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).
 
## [Unreleased] - yyyy-mm-dd
 
### To Be Added

- Coloring of protein sequences based on amino acid properties?
    * If you use Notepad++, I added [a simple tool](https://github.com/molsonkiko/NPP_protein_lexer) that colorizes the amino acids of proteins.
- Easier comparison of motifs and notable features present in one isoform and missing from others?
 
### To Be Changed

- Consider allowing one-line display of protein sequence. Not sure why I'd want to do this.
 
### To Be Fixed

- nothing known

## [0.11.2] - 2022-12-01

1. Make protein list on [index page](https://alignisoforms-production.up.railway.app/) prettier.

## [0.11.1] - 2022-12-01

1. For a protein's page, if no alignment exists, show a button that allows the user to retry the request for a multiple sequence alignment to the EBI computer.

## [0.11.0] - 2022-11-23

1. Add admin form for adding peptides from CSV file.
2. When highlighting all peptides in a protein or alignment, automatically scroll until the first peptide highlighted is in view.

## [0.10.4] - 2022-11-17

1. Tests now work.
2. Better site map.
3. More informative errors when the [get_protein](https://alignisoforms-production.up.railway.app/get_protein/) page fails.

## [0.10.0] - 2022-11-17

1. __First minor release after successful deployment.__
2. Added [about](https://alignisoforms-production.up.railway.app/about) page explaining the motivation of the project.

## [0.9.0] - 2022-11-17

1. Get ready for deployment.
2. Add tests.
3. Remove unnecessary /peptides from the beginning of each address, so the index is the first view you come to, and e.g. `/proteins/P56856` takes you to a protein rather than `/peptides/proteins/P56856`.

## [0.8.0] - 2022-11-14

1. Add data validation for accession number on `peptides\get_protein` form. Now must conform to the regular expression found [here](https://www.uniprot.org/help/accession_numbers).

## [0.7.0] - 2022-11-14

1. Fix bugs in database.
2. Further enhance proteins and alignments views. Can adjust number of amino acids per line.
3. Allow highlighting of all peptides simulaneously (button in alignment and protein form).
4. Add some basic tests.

## [0.6.0] - 2022-11-13

1. Greatly improve appearance of alignments page, allow jumping to mass spec peptides within the alignment while still making the alignment pretty easily readable. Still have some kinks to work out.
2. Protein is now all in one line.
3. Can reorder proteins in the index view. See site_map.html.

## [0.5.0] - 2022-11-13

1. Move data over to a PostgreSQL database.
2. Add site map.
3. Add ability to jump to and highlight a selected mass spec peptide in a sequence on the protein page.

## [0.4.0] - 2022-11-12

1. Improved the website's protein view, adding location data for where each mass spec peptide is found in its protein.
2. Added a simple REST API that takes the UniProt accession number of a protein and returns the isoform IDs, sequence, multiple sequence alignment with isoforms, and mass spec peptide locations in its sequence. See [protein_json_example.json](/protein_json_example.json). The schema for this JSON is at [protein_json_schema.json](/website/peptides/static/peptides/protein_json_schema.json).
3. Add the ability to download the alignment as a `clustal_num` text file from the alignment view.

## [0.3.0] - 2022-11-12

1. Added a very simple, not-yet-production-ready website for including the peptides that Les found in his research and getting the isoforms and a multiple sequence alignment for any protein.

## [0.2.0] - 2022-11-04

1. Changed `align_isoforms(accession_num)` so that it always outputs a CLUSTALO format file for any protein
    with 2+ isoforms, so that proteins with 2 isoforms have the same output format as proteins with 3+ isoforms.

## [0.1.0] - 2022-11-04

### Added

1. Basic API for getting data from UniProt for all isoforms of a protein:
    * `get_all_prots(accession_num)` returns the [UniProt API](https://rest.uniprot.org/docs/#/uniprotkb/searchCursor)
        JSON for all isoforms of the protein with that UniProt accession number.
    * `align_isoforms(accession_num)` gets the sequences of all isoforms of a protein from UniProt,
        and then does a multiple-sequence alignment of those isoforms (if there is more than one).
        Outputs a CLUSTALO format string (e.g., `ampk_gamma.clustal_num`) for multiple-sequence alignment of 3+ isoforms,
        and a simple pairwise sequence alignment string (e.g., `claudin_alignment.txt`) for 2 isoforms.