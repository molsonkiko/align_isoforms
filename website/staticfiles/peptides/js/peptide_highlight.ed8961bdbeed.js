function highlight_locs(loc) {
    // make the most recently highlighted peptide red and bold.
    // undo the highlighting of all other peptides
    var peptides = document.getElementsByClassName('peptide');
    for (var ii = 0; ii < peptides.length; ii++) {
        peptides[ii].style.color = 'black';
        peptides[ii].style.fontWeight = 'normal';
    }
    var new_loc_spans = document.getElementsByClassName(loc);
    for (var ii = 0; ii < new_loc_spans.length; ii++) {
        new_loc_spans[ii].style.color = 'red';
        new_loc_spans[ii].style.fontWeight = 'bold';
    }
}

function highlight_all_peptides() {
    // make all peptides red and bold
    var peptides = document.getElementsByClassName('peptide');
    for (var ii = 0; ii < peptides.length; ii++) {
        peptides[ii].style.color = 'red';
        peptides[ii].style.fontWeight = 'bold';
    }
}

// TODO: Consider reimplementing sequence_chunker in JavaScript 
//     so that all of the re-chunking is done server-side.
//     If I do it that way, I can have the width of proteins
//     and alignments react when the user resizes their browser window.