import os
from natsort import natsorted
from Bio import SeqIO


'''By James C. Hu
This script will:
1) Search through each hadv genotype directory to determine if there are more than 1 gene fasta.
2) Merge the fasta files into a single fasta for MSA.
'''

home_path = os.getcwd()
penton_dir_path = f'{home_path}/clean_gb_files/penton_fastas'
hexon_dir_path = f'{home_path}/clean_gb_files/hexon_fastas'
fiber_dir_path = f'{home_path}/clean_gb_files/fiber_fastas'

penton_genotype_directories = [directory for directory in os.listdir(penton_dir_path)
                               if directory.startswith('penton')]
penton_genotype_directories = natsorted(penton_genotype_directories)
hexon_genotpye_directories = [directory for directory in os.listdir(hexon_dir_path)
                              if directory.startswith('hexon')]
hexon_genotpye_directories = natsorted(hexon_genotpye_directories)
fiber_genotpye_directories = [directory for directory in os.listdir(fiber_dir_path)
                              if directory.startswith('fiber')]
fiber_genotpye_directories = natsorted(fiber_genotpye_directories)


penton_multifile_dirs = [directory for directory in penton_genotype_directories
                         if len(os.listdir(os.path.join(penton_dir_path, directory))) > 1]
hexon_multifile_dirs = [directory for directory in hexon_genotpye_directories
                        if len(os.listdir(os.path.join(hexon_dir_path, directory))) > 1]
fiber_multifile_dirs = [directory for directory in fiber_genotpye_directories
                        if len(os.listdir(os.path.join(fiber_dir_path, directory))) > 1]

os.makedirs('clean_gb_files/combined_fastas/combined_penton_fastas_by_type')
os.makedirs('clean_gb_files/combined_fastas/combined_hexon_fastas_by_type')
os.makedirs('clean_gb_files/combined_fastas/combined_fiber_fastas_by_type')

for directory in penton_multifile_dirs:
    geneotype = directory.split('_')[1]
    with open(f'clean_gb_files/combined_fastas/combined_penton_fastas_by_type/combined_penton_genotype_{geneotype}_sequences.fasta', 'w') as f:
        f.write('')
    fastas = [file for file in os.listdir(f'{penton_dir_path}/{directory}')
              if file.endswith('penton.fasta')]

    for fasta in fastas:
        for seq_record in SeqIO.parse(f'{penton_dir_path}/{directory}/{fasta}', 'fasta'):
            with open(f'clean_gb_files/combined_fastas/combined_penton_fastas_by_type/combined_penton_genotype_{geneotype}_sequences.fasta', 'a') as f:
                f.write(f'>{seq_record.id}\n')
                f.write(f'{seq_record.seq}\n')


for directory in hexon_multifile_dirs:
    geneotype = directory.split('_')[1]
    with open(f'clean_gb_files/combined_fastas/combined_hexon_fastas_by_type/combined_hexon_genotype_{geneotype}_sequences.fasta', 'w') as f:
        f.write('')
    fastas = [file for file in os.listdir(f'{hexon_dir_path}/{directory}')
              if file.endswith('hexon.fasta')]

    for fasta in fastas:
        for seq_record in SeqIO.parse(f'{hexon_dir_path}/{directory}/{fasta}', 'fasta'):
            with open(f'clean_gb_files/combined_fastas/combined_hexon_fastas_by_type/combined_hexon_genotype_{geneotype}_sequences.fasta', 'a') as f:
                f.write(f'>{seq_record.id}\n')
                f.write(f'{seq_record.seq}\n')


for directory in fiber_multifile_dirs:
    geneotype = directory.split('_')[1]
    with open(f'clean_gb_files/combined_fastas/combined_fiber_fastas_by_type/combined_fiber_genotype_{geneotype}_sequences.fasta', 'w') as f:
        f.write('')
    fastas = [file for file in os.listdir(f'{fiber_dir_path}/{directory}')
              if file.endswith('fiber.fasta')]

    for fasta in fastas:
        for seq_record in SeqIO.parse(f'{fiber_dir_path}/{directory}/{fasta}', 'fasta'):
            with open(f'clean_gb_files/combined_fastas/combined_fiber_fastas_by_type/combined_fiber_genotype_{geneotype}_sequences.fasta', 'a') as f:
                f.write(f'>{seq_record.id}\n')
                f.write(f'{seq_record.seq}\n')
