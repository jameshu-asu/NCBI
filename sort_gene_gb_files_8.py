import os
import shutil
from natsort import natsorted
from Bio import SeqIO

'''By James C. Hu
This script will:
1) Consolidate all fastas of the same gene into a single fasta file.
2) Generate directoreis for all identified hadv genotypes
3) Move all the individual fasta files into their respective directories.
'''


def grab_genotype(fasta_list: list) -> list:
    '''This function will:
    Create a list of hadv genotypes using the fasta file names.
    '''
    genotype_list = []
    for fasta in fasta_list:
        if fasta.startswith(('AC', 'NC')):
            genotype = fasta.split('_')[4]
            genotype_list.append(genotype)
        else:
            genotype = fasta.split('_')[3]
            genotype_list.append(genotype)
    return genotype_list


def sort_fastas(fasta_list: list, directory: str) -> None:
    '''This function will:
    Move specific hadv files into their respective directories.

    directory: case sensitive, only takes 'penton_fastas', 'hexon_fastas', or 'fiber_fastas'
    '''
    for fasta in fasta_list:
        if fasta.startswith(('AC', 'NC')):
            genotype = fasta.split('_')[4]
        else:
            genotype = fasta.split('_')[3]
        shutil.move(f'clean_gb_files/{directory}/{fasta}',
                    f'clean_gb_files/{directory}/{directory[:-7]}_{genotype}')
    return None


os.makedirs('clean_gb_files/combined_fastas')

penton_fasta_list = [file for file in os.listdir('clean_gb_files/penton_fastas')
                     if file.endswith('penton.fasta')]
penton_fasta_list = natsorted(list(set(penton_fasta_list)))

hexon_fasta_list = [file for file in os.listdir('clean_gb_files/hexon_fastas')
                    if file.endswith('hexon.fasta')]
hexon_fasta_list = natsorted(list(set(hexon_fasta_list)))

fiber_fasta_list = [file for file in os.listdir('clean_gb_files/fiber_fastas')
                    if file.endswith('fiber.fasta')]
fiber_fasta_list = natsorted(list(set(fiber_fasta_list)))


with open(f'clean_gb_files/combined_fastas/combined_all_penton_sequences.fasta', 'w') as f:
    pass

with open(f'clean_gb_files/combined_fastas/combined_all_hexon_sequences.fasta', 'w') as f:
    pass

with open(f'clean_gb_files/combined_fastas/combined_all_fiber_sequences.fasta', 'w') as f:
    pass


for file in penton_fasta_list:
    for seq_record in SeqIO.parse(f'clean_gb_files/penton_fastas/{file}', 'fasta'):
        with open(f'clean_gb_files/combined_fastas/combined_all_penton_sequences.fasta', 'a') as f:
            f.write(f'>{seq_record.id}\n')
            f.write(f'{seq_record.seq}\n')

for file in hexon_fasta_list:
    for seq_record in SeqIO.parse(f'clean_gb_files/hexon_fastas/{file}', 'fasta'):
        with open(f'clean_gb_files/combined_fastas/combined_all_hexon_sequences.fasta', 'a') as f:
            f.write(f'>{seq_record.id}\n')
            f.write(f'{seq_record.seq}\n')

for file in fiber_fasta_list:
    for seq_record in SeqIO.parse(f'clean_gb_files/fiber_fastas/{file}', 'fasta'):
        with open(f'clean_gb_files/combined_fastas/combined_all_fiber_sequences.fasta', 'a') as f:
            f.write(f'>{seq_record.id}\n')
            f.write(f'{seq_record.seq}\n')

penton_genotype_list = natsorted(list(set(grab_genotype(penton_fasta_list))))
hexon_genotype_list = natsorted(list(set(grab_genotype(hexon_fasta_list))))
fiber_genotype_list = natsorted(list(set(grab_genotype(fiber_fasta_list))))

all_genotype_list = natsorted(list(set(penton_genotype_list
                                       + hexon_genotype_list
                                       + fiber_genotype_list)))


for genotype in all_genotype_list:
    try:
        os.mkdir(f'clean_gb_files/penton_fastas/penton_{genotype}')
    except Exception as e:
        print(e)
        pass
    try:
        os.mkdir(f'clean_gb_files/hexon_fastas/hexon_{genotype}')
    except Exception as e:
        print(e)
        pass
    try:
        os.mkdir(f'clean_gb_files/fiber_fastas/fiber_{genotype}')
    except Exception as e:
        print(e)
        pass

sort_fastas(penton_fasta_list, 'penton_fastas')
sort_fastas(hexon_fasta_list, 'hexon_fastas')
sort_fastas(fiber_fasta_list, 'fiber_fastas')
