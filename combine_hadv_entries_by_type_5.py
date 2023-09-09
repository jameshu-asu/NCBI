import pandas as pd
import os
from Bio import SeqIO

'''By James C. Hu
This script will:
Use the final post qc accession list to extract the desired entries from the original gb file.
'''

df = pd.read_csv('qc_logs/final_accession_list.tsv', sep='\t')

log_file = 'clean_gb_files/logs/directory_setup_log.txt'
gb_file = 'initial_download/genbank_files/deduped_combined_NCBI_hadv_entries.gb'

dict1 = {}

os.makedirs(f'clean_gb_files/HAdV_complete_genome')
os.makedirs(f'clean_gb_files/logs')

with open(log_file, 'w') as f:
    f.write(f'Genotype\tRecord_count\n')

for group in df.groupby('Manual_type_check'):
    genotype = group[0]
    accession_list = group[1]['Accession'].to_list()
    with open(log_file, 'a') as f:
        f.write(f'{genotype}\t{len(accession_list)}\n')
    dict1[genotype] = accession_list

for genotype in dict1.keys():
    with open(f'clean_gb_files/HAdV_complete_genome/HAdV_{genotype}.gb', 'w') as f:
        pass

for entry in SeqIO.parse(gb_file, 'gb'):
    accession = entry.name
    seq = entry.seq
    for genotype, accession_list in dict1.items():
        if accession in accession_list:
            with open(f'clean_gb_files/HAdV_complete_genome/HAdV_{genotype}.gb', 'a') as f:
                SeqIO.write(entry, f, 'gb')
