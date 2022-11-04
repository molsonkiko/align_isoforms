# Change Log
All [notable changes](#420---2022-10-29) to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).
 
## [Unreleased] - yyyy-mm-dd
 
### To Be Added

- Coloring of protein sequences based on amino acid properties?
    * If you use Notepad++, I added [a simple tool](https://github.com/molsonkiko/NPP_protein_lexer) that colorizes the amino acids of proteins.
- Easier comparison of motifs and notable features present in one isoform and missing from others?
- License for repository (what kind to use? MIT? Apache?)
 
### To Be Changed

- nothing known
 
### To Be Fixed

- nothing known

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