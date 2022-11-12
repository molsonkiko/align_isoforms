'''Parse the file containing all the unique peptides
found for each variant
'''
import json
import pandas as pd

df = pd.read_csv('unique_peptides.csv')

def uniques_per_acc_num(df: pd.DataFrame) -> dict:
    out = {}
    for _, row in df.iterrows():
        acc_num = row.accession_number
        pep = row.peptide
        out.setdefault(acc_num, set())
        other_peptides = out[acc_num]
        has_super_peptide = False
        for other_pep in list(other_peptides):
            if pep in other_pep:
                has_super_peptide = True
            if other_pep in pep:
                other_peptides.remove(other_pep)
        if not has_super_peptide:
            other_peptides.add(pep)
    return {k: list(v) for k, v in out.items()}

if __name__ == '__main__':
    uniques = uniques_per_acc_num(df)
    with open('unique_peptides_per_acc_num.json', 'w') as f:
        json.dump(uniques, f, indent=4)