let last_highlighted_locs = "0";

function highlight_locs(loc) {
    // make the most recently highlighted peptide red and bold.
    // undo the highlighting of the previously highlighted peptide
    // (if any)
    var last_loc_spans = document.getElementsByClassName(last_highlighted_locs);
    if (last_loc_spans !== undefined) {
        for (var ii = 0; ii < last_loc_spans.length; ii++) {
            last_loc_spans[ii].style.color = 'black';
            last_loc_spans[ii].style.fontWeight = 'normal';
        }
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
    for (let ii = 0; ii < peptides.length; ii++) {
        peptides[ii].style.color = 'red';
        peptides[ii].style.fontWeight = 'bold';
    }
}