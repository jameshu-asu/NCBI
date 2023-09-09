import pandas as pd
from math import nan

'''
This script will:
1) Attempt to assing hadv genotype by matching data in the 'Source' field, with data from the 'Organism' field.
2) Assign a temporary value for 'Source', 'Organism' missmatch.
3) Attempt to assign type using data from the 'Strain' field using regex.
4) Attempt to extract hadv type data from 'Title' field by removing all alphabetical characters using regex
5) Manually assing types for samples that failed to be assigned during steps 1-4
    - These samples may have had type data logged elseware in the genbank file.
    - These samples may have had conflicting information in the above fields as a result of the regex based character removal process.
6) Manually remove genbank entries with ambiguous hadv type data.
    - These entries had either conflicting genotype information or genotype data was missing all together.
'''

# pd options to adjust terminal display
pd.set_option('display.max_columns', 15)
pd.set_option('display.width', 1000)

df = pd.read_csv('qc_logs/missmatched_output_1.tsv', sep='\t').copy()
df = df.fillna(0)
df.insert(7, 'Manual_type_check', nan)
df[['Source', 'Organism', 'Serotype']] = df[[
    'Source', 'Organism', 'Serotype']].astype(int)

temp_df1 = pd.DataFrame()
temp_df2 = pd.DataFrame()
temp_df3 = pd.DataFrame()
temp_df4 = pd.DataFrame()

#____Manual_type_check___#
# Assigning type as Source if Source == Organism
temp_df1['Accession'] = df['Accession']
temp_df1['Manual_type_check'] = df['Source'].loc[df['Source'] == df['Organism']]
temp_df1 = temp_df1.set_index('Accession')

# Assigning type as Serotype if data is availabe, and Source and Organism both == 0
temp_df2['Accession'] = df['Accession']
temp_df2['Manual_type_check'] = df['Serotype'].loc[(
    df['Source'] == 0) & (df['Organism'] == 0)]
temp_df2 = temp_df2.set_index('Accession')

# Assigning type based on Strain data using regex
temp_df3['Accession'] = df['Accession']
temp_df3['Manual_type_check'] = df['Strain'].str.extract(r'(\d+)\[[P]\d+[H]\d+[F]\d+]')
temp_df3 = temp_df3.fillna(nan)
temp_df3 = temp_df3.set_index('Accession')

temp_df1.replace(0, nan, inplace=True)
temp_df2.replace(0, nan, inplace=True)

temp_df1.update(temp_df2)
temp_df1.update(temp_df3)

df = df.set_index('Accession')
df.update(temp_df1)

# Extracting type data from titles
temp_df4['Manual_type_check'] = temp_df1.loc[temp_df1['Manual_type_check'].isna()]
# temp_df4_accession_list: list of accessions with nan in Manual_type_check field.
temp_df4_accession_list = temp_df4.index.to_list()
temp_df4 = df.loc[[i for i in temp_df4_accession_list]].copy()
temp_df4 = temp_df4[['Manual_type_check', 'Titles']]
temp_df4['Manual_type_check'] = temp_df4['Titles'].str.extract(r'(\d+)')
df.update(temp_df4)
df.to_csv('qc_logs/missmatched_output_2.tsv', sep='\t')

df_in1 = pd.read_csv('qc_logs/qc5_output.tsv', sep='\t').copy()
df_in1.insert(7, 'Manual_type_check', nan)
df_in1['Manual_type_check'] = df_in1['Source']
df_in2 = pd.read_csv('qc_logs/missmatched_output_2.tsv',
                     sep='\t')[['Accession', 'Manual_type_check']].copy()
df_in1 = df_in1.set_index('Accession')
df_in2 = df_in2.set_index('Accession')
df_in1.update(df_in2)

# manually setting type data for edge case annotations.
df_in1.at['MN531562', 'Manual_type_check'] = '7'
df_in1.at['KF268122', 'Manual_type_check'] = '37'
df_in1.at['MK501726', 'Manual_type_check'] = '53'

# removing post qc accessions with ambiguous type data.
df_in1 = df_in1[df_in1['Manual_type_check'].notna()]
df_in1 = df_in1[df_in1['Manual_type_check'].astype(int) < 113]
df_in1.to_csv('qc_logs/missmatched_output_3.tsv', sep='\t')

# generating final accession/hadv type list
final_accession_df = df_in1[['Manual_type_check']].astype(int)
final_accession_df = final_accession_df.sort_values('Manual_type_check')
final_accession_df.to_csv('qc_logs/final_accession_list.tsv', sep='\t')
