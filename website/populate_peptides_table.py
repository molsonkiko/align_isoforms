import pandas as pd
import sqlite3

df = (pd.read_csv('peptides/static/peptides/unique_peptides_per_acc_num.csv')
    .rename({'acc_num': 'prot', 'peps': 'peptide'}, axis = 1))

conn = sqlite3.connect('db.sqlite3')

df.to_sql('peptides_peptide', conn, if_exists='append', index=False)