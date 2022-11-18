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

There is now [a website](https://alignisoforms-production.up.railway.app/) that can do all of these things and has a handy tool for visualizing the locations of mass spec peptides of interest.

Current version: [0.10.4](/CHANGELOG.md#0104---2022-11-17)

![Example of alignment with highlighted peptides](/website/alignment%20page%20all_peptides_highlighted.PNG)