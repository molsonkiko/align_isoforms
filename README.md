# align_isoforms
Working in collaboration with biology professor [Les Timpe](https://github.com/lestimpe) to develop a tool for biologists to easily compare the primary sequences of different isoforms of proteins.


Basic use:
* `get_all_prots(accession_num)` returns the [UniProt API](https://rest.uniprot.org/docs/#/uniprotkb/searchCursor)
    JSON for all isoforms of the protein with that UniProt accession number.
* `align_isoforms(accession_num)` gets the sequences of all isoforms of a protein from UniProt,
    and then does a multiple-sequence alignment of those isoforms (if there is more than one).
    Outputs a CLUSTALO format string (e.g., `ampk_gamma.clustal_num`, `claudin18.clustal_num`)
    for multiple-sequence alignment of any protein with 2+ isoforms.
    Outputs the sequence as a string for any proteins with only one isoform.

I am also working on a website that can do all these things.

Current version: [0.3.0](/CHANGELOG.md#030---2022-11-12)