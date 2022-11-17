import pandas as pd
from sqlalchemy import create_engine

df = (pd.read_csv('peptides/static/peptides/unique_peptides_per_acc_num.csv')
    .rename({'acc_num': 'prot', 'peps': 'peptide'}, axis = 1))
df['location'] = -1

username = 'postgres'
password = input('Enter the password for the postgres database: ')
host = input('Enter the host: ')
port = input('Enter the port: ')
dbname = 'railway'
conn_str = f'postgresql+psycopg2://{username}:{password}@{host}:{port}/{dbname}'

engine = create_engine(conn_str)
with engine.connect() as conn:
    df.to_sql('peptides_peptide', conn, if_exists='append', index=False)
conn.close()
engine.dispose()