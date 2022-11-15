import pandas as pd
from sqlalchemy import create_engine

df = (pd.read_csv('peptides/static/peptides/unique_peptides_per_acc_num.csv')
    .rename({'acc_num': 'prot', 'peps': 'peptide'}, axis = 1))
df['location'] = -1

port = '5433'
username = 'postgres'
dbname = 'splice_variants_website'
with open(r'..\..\svs\postgres_password.txt') as f:
    password = f.read()
conn_str = f'postgresql+psycopg2://postgres:{password}@localhost:{port}/{dbname}'

engine = create_engine(conn_str)
with engine.connect() as conn:
    df.to_sql('peptides_peptide', conn, if_exists='append', index=False)
conn.close()
engine.dispose()