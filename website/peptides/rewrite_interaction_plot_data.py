import re
import os
import pathlib

if __name__ == '__main__':
    FDIR = pathlib.Path(__file__).parent/'static'/'peptides'/'isoform abundance cancer vs not'
    for fname in os.listdir(FDIR):
        if '.' in fname:
            continue
        old_fname = FDIR / fname
        with open(old_fname) as f:
            lines = re.split('\r?\n', f.read())
            lines[0] = lines[0].replace('"', '').replace(' ', ',')
            for ii in range(1, len(lines)):
                lines[ii] = re.sub('^"\d+" ', '', lines[ii]).replace(' ', ',')
        new_fname = FDIR / (fname + '.csv')
        with new_fname.open('w') as f:
            f.write('\n'.join(lines))
        old_fname.unlink()