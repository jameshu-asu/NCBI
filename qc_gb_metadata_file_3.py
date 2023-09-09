import pandas as pd
from math import nan
import os

'''By James C. Hu
This script will:
1) Filter gb entries by match criteria
    str_match1 = 'Human|Homo|human|homo'
    str_match2 = 'Adenovirus|Mastadenovirus|adenovirus|mastadenovirus'https://github.com/jameshu-asu/NCBI/blob/master/qc_gb_metadata_file_3.py
2) Filter gb entries by genome length
    min_genomoe_length = 30_000
3) Deduplicate gb entries
4) Extract Strain data using regex
5) Strip other categoris of letter to try and extract HAdV genotype number using regex
6) Generate a log file for each process.
'''

# pd options to adjust terminal display
pd.set_option('display.max_column', 15)
pd.set_option('display.width', 1000)


def remove_chars(df: pd.DataFrame, column_name: str) -> list:
    '''
    This function will replace all string character values within a series object with a spaces using regex.
    '''
    df[column_name] = df[column_name].str.replace('\D', '', regex=True)
    source_list = [i for i in df[column_name].to_list()]
    for idx, item in enumerate(source_list):
        if item == '':
            source_list[idx] = nan
    return source_list


# make sure this value makes sense.
min_genomoe_length = 30_000

os.mkdir('qc_logs')

df = pd.read_csv('initial_download/download_logs/hadv_genbank_metadata.tsv', sep='\t').copy()

# string filters, selecting FOR these strings
str_match1 = 'Human|Homo|human|homo'
str_match2 = 'Adenovirus|Mastadenovirus|adenovirus|mastadenovirus'

# filtering by match1
df_qc1 = df.loc[df['Source'].str.contains(str_match1)].copy()

# filtering by match2
df_qc2 = df_qc1.loc[df['Source'].str.contains(str_match2)].copy()

# filtering by genome length
df_qc3 = df_qc2.loc[df['Ref_genome_length'] > min_genomoe_length].copy()

# deduplicating accession ids
df_qc4 = df_qc3.drop_duplicates(subset='Accession', keep='first').copy()

# extracting penton, hexon and fiber alleles
df_qc4['P_H_F_allele'] = df_qc4['Strain'].str.extract(
    r'([P]\d+[H]\d+[F]\d+)')

df_qc4['Penton_allele'] = df_qc4['P_H_F_allele'].str.extract(r'P(\d+)')
df_qc4['Hexon_allele'] = df_qc4['P_H_F_allele'].str.extract(r'H(\d+)')
df_qc4['Fiber_allele'] = df_qc4['P_H_F_allele'].str.extract(r'F(\d+)')

# extracting numeric type data from 'Source', 'Organism', and 'Serotype' columns.
df_qc5 = df_qc4.copy()
df_qc5['Source'] = remove_chars(df_qc5, 'Source')
df_qc5['Organism'] = remove_chars(df_qc5, 'Organism')
df_qc5['Serotype'] = remove_chars(df_qc5, 'Serotype')
df_qc5['Serotype'].fillna(nan, inplace=True)

# extracting accessions that were incorrectly grouped during the initial download.
df_qc5_temp1 = df_qc5.loc[df_qc5['Source'].isnull()].copy()
df_qc5_temp2 = df_qc5.loc[df_qc5['Source'].notna()].copy()
df_qc5_temp2['Source'] = df_qc5_temp2['Source'].astype(int)
df_qc5_temp2 = df_qc5_temp2[df_qc5_temp2['Source'] > 113]

df_missmatched = pd.concat([df_qc5_temp1, df_qc5_temp2])
df_missmatched = df_missmatched.sort_values('Accession')

df_gene_alleles = df_qc5[['Accession', 'P_H_F_allele', 'Penton_allele',
                          'Hexon_allele', 'Fiber_allele']]
df_gene_alleles = df_gene_alleles.dropna(axis=0)

# generate output logs
df_qc1.to_csv('qc_logs/qc1_output.tsv', sep='\t', index=None)
df_qc2.to_csv('qc_logs/qc2_output.tsv', sep='\t', index=None)
df_qc3.to_csv('qc_logs/qc3_output.tsv', sep='\t', index=None)
df_qc4.to_csv('qc_logs/qc4_output.tsv', sep='\t', index=None)
df_qc5.to_csv('qc_logs/qc5_output.tsv', sep='\t', index=None)
df_missmatched.to_csv('qc_logs/missmatched_output_1.tsv', sep='\t', index=None)
df_gene_alleles.to_csv('qc_logs/hadv_allele_data.tsv', sep='\t', index=None)
