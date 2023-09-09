'''
Dependencies:
Pandas 2.0+
Biopython 1.80
natsort 8.3.1
'''

import pandas as pd
import Bio
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio import SeqIO
from itertools import zip_longest
from natsort import natsorted
import os
from datetime import datetime

'''By James C. Hu
This script will:
1) Attempt to extract all relevant product names from each CDS section.
2) Filter the list of product names.
3) Use the product names to extract gene loci by nt position
4) Extract and export the specific gene sequence using the locis from 3.
'''

# adjust pandas terminal output
pd.options.display.max_columns = 10
pd.options.display.width = 1000


def seq_start(string: str) -> int:
    '''
    This function strips the input sequence location string and returns the first number (starting nucleotide).
    Locations that start with 'join' (genes that may have an upstream promoter region) are skipped.
    '''
    if str(string).startswith('join'):
        return None
    else:
        return int(str(string).strip('[]').strip('(+)').strip('(-)').strip(']').split(':')[0])


def seq_end(string: str) -> int:
    '''
    This function strips the input sequence location string and returns the last number (ending nucleotide).
    Locations that start with 'join' (genes that may have an upstream promoter region) are skipped.
    '''
    if str(string).startswith('join'):
        return None
    else:
        return int(str(string).strip('[]').strip('(+)').strip('(-)').strip(']').split(':')[1])


def pull_SimpleLocation_seq(seq_start: int, seq_end: int) -> str:
    '''
    This function extracts sequences from SeqIO.Seq objects based on nucleotide start and stop postitions.
    '''
    try:
        return SeqFeature(SimpleLocation(seq_start, seq_end), type='CDS').extract(seq)
    except Exception:
        return None


def feature_selection(record: Bio.SeqRecord.SeqRecord) -> list:
    '''
    This function searches the gb file as a Bio.SeqRecord.SeqRecord object and returns all product names for all features that are NOT blank as a deduplicated flat list
    '''
    feature_list = []
    for feature in record.features:
        try:
            feature_list.append([feature for feature in feature.qualifiers.get(
                'product') if feature != 'NoneType'])
        except Exception:
            pass
    unique_flat_feature_list = list(set([j for i in feature_list for j in i]))
    return unique_flat_feature_list


def grab_HAdV_qualifers(file_name: str) -> None:
    '''
    This function writes the gene id lists to a text file.
    Must define relevant lists before calling function
    '''
    with open(file_name, 'w') as f:
        f.write(f'Penton_qualifiers\tHexon_qualifiers\tFiber_qualifiers\n')
        for penton, hexon, fiber in zip_longest(penton_id_list, hexon_id_list, fiber_id_list):
            f.write(f'{penton}\t{hexon}\t{fiber}\n')
    return None


os.mkdir(f'clean_gb_files/penton_fastas')
os.mkdir(f'clean_gb_files/hexon_fastas')
os.mkdir(f'clean_gb_files/fiber_fastas')

penton_id_list = []
hexon_id_list = []
fiber_id_list = []

# Search current directory for files. Can specify path by passing path into os.listdir()
gb_in_file = 'initial_download/genbank_files/deduped_combined_NCBI_hadv_entries.gb'
# sorts files where character space for digits are not the same. Requires natsort python library


infile_qualifier_ids = 'qualifier_ids.csv'
out_file = 'clean_gb_files/logs/HAdV_gene_location_metadata.txt'

df_qualifiers = pd.read_csv(infile_qualifier_ids).copy()


for seq_record in SeqIO.parse(gb_in_file, 'gb'):
    feature_list = feature_selection(seq_record)
    penton_feature_ids = [qualifier for qualifier in df_qualifiers['Penton_qualifiers'].to_list()
                          if str(qualifier) != 'nan']
    hexon_feature_ids = [qualifier for qualifier in df_qualifiers['Hexon_qualifiers'].to_list()
                         if str(qualifier) != 'nan']
    fiber_feature_ids = [qualifier for qualifier in df_qualifiers['Fiber_qualifiers'].to_list()
                         if str(qualifier) != 'nan']

    for feature in feature_list:
        if any(i in feature for i in penton_feature_ids):
            penton_id_list.append(feature)
    for feature in feature_list:
        if any(i in feature for i in hexon_feature_ids):
            hexon_id_list.append(feature)
    for feature in feature_list:
        if any(i in feature for i in fiber_feature_ids):
            fiber_id_list.append(feature)


penton_id_list = list(set(penton_id_list))
hexon_id_list = list(set(hexon_id_list))
fiber_id_list = list(set(fiber_id_list))

# Remove lists contain string identifiers used to remove unwanted qualifiers from
penton_remove_ids = ['pIIIa', 'IIIa' 'pVII', 'pV', 'pX', 'X', 'V', 'pmu', 'core', 'precursor', 'Mu',
                     'associated', 'hexon-associated', 'tegument', 'UL29', 'UL2', 'glycoprotein',
                     'Pre-hexon-linking', 'pre-hexon-linking', 'L1', 'protein IIIa']

hexon_remove_ids = ['IIIa', 'IX', 'pIX', 'pVIII', 'VIII', '100', '100kDa', 'L2', 'glycoprotein',
                    'associated', 'hexon-associated', 'proteinase', 'protease', 'UL30A', 'hexon-assembly',
                    'Protease', 'precursor', 'hexon-interlacing', 'Hexon-interlacing', 'UL33', 'UL36', 'UL38',
                    'penton', 'UL30', 'tegument', 'UL34', 'pIII', 'UL32', 'E1A', 'UL31', 'protein VI', 'pVI', 'VA RNA II',
                    'U exon', 'U-exon', 'U Exon', 'U-Exon']

fiber_remove_ids = ['proteincomplement', 'agnoprotein', 'RL5A', 'tegument', 'IVA2', 'pIVA2', 'maturation',
                    'encapsidation', 'IVa2', 'UL5']

updated_penton_id_list = []
updated_hexon_id_list = []
updated_fiber_id_list = []

# removes unwanted qualifiers from gene id list.
for penton_id in penton_id_list:
    if all(substring not in penton_id for substring in penton_remove_ids):
        updated_penton_id_list.append(penton_id)
# print(len(updated_penton_id_list))
# print(updated_penton_id_list)

# print(len(hexon_id_list))
# print(hexon_id_list)
for hexon_id in hexon_id_list:
    if all(substring not in hexon_id for substring in hexon_remove_ids):
        updated_hexon_id_list.append(hexon_id)
# print(len(updated_hexon_id_list))
# print(updated_hexon_id_list)

# print(len(fiber_id_list))
# print(fiber_id_list)
for fiber_id in fiber_id_list:
    if all(substring not in fiber_id for substring in fiber_remove_ids):
        updated_fiber_id_list.append(fiber_id)
# print(len(updated_fiber_id_list))
# print(updated_fiber_id_list)

penton_qalifier_list = []
hexon_qalifier_list = []
fiber_qalifier_list = []

for penton_id in updated_penton_id_list:
    penton_id_list = [penton_id]
    penton_qalifier_list.append(penton_id_list)

for hexon_id in updated_hexon_id_list:
    hexon_id_list = [hexon_id]
    hexon_qalifier_list.append(hexon_id_list)

for fiber_id in updated_fiber_id_list:
    fiber_id_list = [fiber_id]
    fiber_qalifier_list.append(fiber_id_list)

# only need to unhash if this specific output is needed again, also unhash lines 240-241.
with open(out_file, 'w') as f:
    f.write('Accession_ID\tVersion_ID\tGB_submission_date\tGenome_size\tOrganism\tTaxonomy\tAuthors\tTitle\tJournal\tPenton_product\tHexon_product\tFiber_product\tPenton_location\tHexon_location\tFiber_location\n')

clean_gb_file_list = [file for file in os.listdir('clean_gb_files/HAdV_complete_genome')
                      if file.startswith('HAdV')]
clean_gb_file_list = natsorted(clean_gb_file_list)

for file in clean_gb_file_list:
    genotype = file.split('_')[1].split('.')[0]
    # print(genotype)
    for record in SeqIO.parse(f'clean_gb_files/HAdV_complete_genome/{file}', 'gb'):
        accession = record.id[:-2]
        version = record.id
        # print(accession)
        date1 = datetime.strptime(record.annotations['date'], '%d-%b-%Y').strftime('%Y-%m-%d')
        size = len(record.seq)
        organism = record.annotations['organism']
        taxonomy = record.annotations['taxonomy']
        reference_title = record.annotations['references'][0].title
        reference_authors = record.annotations['references'][0].authors
        reference_journal = record.annotations['references'][0].journal
        penton_location = ''
        hexon_location = ''
        fiber_location = ''
        DNA_pol_location = ''
        penton_product = ''
        hexon_product = ''
        fiber_product = ''

        seq = record.seq

        for feature in record.features:
            # Notes: qualifiers are returned as lists
            if feature.qualifiers.get('product') in penton_qalifier_list:
                penton_product = feature.qualifiers.get('product')
                penton_location = feature.location
                penton_start = seq_start(penton_location)
                penton_end = seq_end(penton_location)
                penton_seq = pull_SimpleLocation_seq(penton_start, penton_end)

                # fasta writing steps can be commented out if they have already been generated.
                # Note: Files will be written over if uncommented.
                with open(f'clean_gb_files/penton_fastas/{accession}_HAdV_type_{genotype}_penton.fasta', 'w') as f:
                    f.write(f'>{accession}_HAdV_type_{genotype}_penton\n')
                    f.write(f'{str(penton_seq)}\n')

            elif feature.qualifiers.get('product') in hexon_qalifier_list:
                hexon_product = feature.qualifiers.get('product')
                hexon_location = feature.location
                hexon_start = seq_start(hexon_location)
                hexon_end = seq_end(hexon_location)
                hexon_seq = pull_SimpleLocation_seq(hexon_start, hexon_end)

                with open(f'clean_gb_files/hexon_fastas/{accession}_HAdV_type_{genotype}_hexon.fasta', 'w') as f:
                    f.write(f'>{accession}_HAdV_type_{genotype}_hexon\n')
                    f.write(f'{str(hexon_seq)}\n')

            elif feature.qualifiers.get('product') in fiber_qalifier_list:
                fiber_product = feature.qualifiers.get('product')
                fiber_location = feature.location
                fiber_start = seq_start(fiber_location)
                fiber_end = seq_end(fiber_location)
                fiber_seq = pull_SimpleLocation_seq(fiber_start, fiber_end)

                with open(f'clean_gb_files/fiber_fastas/{accession}_HAdV_type_{genotype}_fiber.fasta', 'w') as f:
                    f.write(f'>{accession}_HAdV_type_{genotype}_fiber\n')
                    f.write(f'{str(fiber_seq)}\n')

            # elif feature.qualifiers.get('product') == ['DNA polymerase']:
            #     DNA_pol_location = feature.location
            # DNA_pol_start = seq_start(DNA_pol_location)
            # DNA_pol_end = seq_end(DNA_pol_location)
            # DNA_pol_seq = pull_seq(DNA_pol_start, DNA_pol_end)
            # print(DNA_pol_location)

        with open(out_file, 'a') as f:
            f.write(f'{accession}\t{version}\t{date1}\t{size}\t{organism}\t{taxonomy}\t{reference_authors}\t{reference_title}\t{reference_journal}\t{penton_product}\t{hexon_product}\t{fiber_product}\t{penton_location}\t{hexon_location}\t{fiber_location}\n')
