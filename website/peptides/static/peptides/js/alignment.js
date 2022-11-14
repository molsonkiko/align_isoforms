let last_highlighted_locs = "0";

function highlight_locs(loc) {
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