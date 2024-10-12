import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from airflow import DAG
from airflow.operators.empty import EmptyOperator
from airflow.operators.python import PythonOperator
from airflow.providers.postgres.hooks.postgres import PostgresHook
from airflow.providers.amazon.aws.hooks.s3 import S3Hook
import pendulum


def extract_data(ti, **kwargs):
    execution_date = kwargs['execution_date']
    ti.xcom_push(key='execution_date', value=execution_date.strftime('%Y-%m-%d'))
    query = """
    SELECT mol_id, name, smiles
    FROM molecules
    LIMIT 100;
    """
    postgres_hook = PostgresHook(postgres_conn_id='postgres_default')
    df = postgres_hook.get_pandas_df(sql=query)
    df.to_csv('/tmp/extracted_data.csv', index=False)


def transform_data():
    df = pd.read_csv('/tmp/extracted_data.csv')

    mol_weights = []
    logP_values = []
    tpsa_values = []
    h_donors = []
    h_acceptors = []
    lipinski_pass = []

    for smiles in df['smiles']:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol_weights.append(Descriptors.MolWt(mol))
            logP_values.append(Descriptors.MolLogP(mol))
            tpsa_values.append(Descriptors.TPSA(mol))
            h_donors.append(Descriptors.NumHDonors(mol))
            h_acceptors.append(Descriptors.NumHAcceptors(mol))

            lipinski = (
                Descriptors.MolWt(mol) < 500 and
                Descriptors.MolLogP(mol) < 5 and
                Descriptors.NumHDonors(mol) <= 5 and
                Descriptors.NumHAcceptors(mol) <= 10
            )
            lipinski_pass.append(lipinski)
        else:
            mol_weights.append(np.nan)
            logP_values.append(np.nan)
            tpsa_values.append(np.nan)
            h_donors.append(np.nan)
            h_acceptors.append(np.nan)
            lipinski_pass.append(False)

    df['Molecular Weight'] = mol_weights
    df['LogP'] = logP_values
    df['TPSA'] = tpsa_values
    df['H Donors'] = h_donors
    df['H Acceptors'] = h_acceptors
    df['Lipinski Pass'] = lipinski_pass

    df.to_excel('/tmp/transformed_data.xlsx', index=False)


def load_data(ti):
    execution_date = ti.xcom_pull(task_ids='extract_data', key='execution_date')
    s3_hook = S3Hook(aws_conn_id='aws_default')
    s3_bucket = 'hw-bucket-plsbeunique'
    s3_key = f'molecule_data/transformed_data_{execution_date}.xlsx'

    s3_hook.load_file(
        filename='/tmp/transformed_data.xlsx',
        key=s3_key,
        bucket_name=s3_bucket,
        replace=True
    )


with DAG(
    dag_id='etl_moleculedata_dag',
    start_date=pendulum.today(),
    schedule_interval='@daily',
    tags=['molecules', 'rdkit', 's3'],
) as dag:

    start_op = EmptyOperator(
        task_id='start'
    )

    extract_data_op = PythonOperator(
        task_id='extract_data',
        python_callable=extract_data,
    )

    transform_data_op = PythonOperator(
        task_id='transform_data',
        python_callable=transform_data,
    )

    load_data_op = PythonOperator(
        task_id='load_data',
        python_callable=load_data,
    )

    finish_op = EmptyOperator(
        task_id='finish'
    )

    start_op >> extract_data_op >> transform_data_op
    transform_data_op >> load_data_op >> finish_op
