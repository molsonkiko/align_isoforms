let last_highlighted_loc = "0";

function highlight_loc(loc) {
    last_loc_span = document.getElementById(last_highlighted_loc);
    if (last_loc_span !== null) {
        document.getElementById(last_highlighted_loc).style.color = 'black';
        document.getElementById(last_highlighted_loc).style.fontWeight = 'normal';
    }
    document.getElementById(loc).style.color = 'red';
    document.getElementById(loc).style.fontWeight = 'bold';
}