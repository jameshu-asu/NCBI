from Bio import SeqIO

'''
This script will:
Compile metadata from the combined human adenovirus (HAdV) gb file.

Metadata logged:
1) Accession ID
2) Reference genome length
3) Source data
4) Organism data
5) Serotype data
6) Acronyms
7) Strain data
8) Notes from the gb file
9) Title of submission if applicable
'''

out_file = 'initial_download/download_logs/hadv_genbank_metadata.tsv'
in_gb_file = 'initial_download/genbank_files/combined_NCBI_hadv_entries.gb'
out_gb_file = 'initial_download/genbank_files/deduped_combined_NCBI_hadv_entries.gb'

# Deduping entries from combined gb file
accession_list = []

with open(out_gb_file, 'w') as f:
    pass

for entry in SeqIO.parse(in_gb_file, 'gb'):
    accession = entry.name
    if accession not in accession_list:
        with open(out_gb_file, 'a') as f:
            SeqIO.write(entry, f, 'gb')
            accession_list.append(accession)
accession_list.clear()

with open(out_file, 'w') as f:
    f.write('Accession\tRef_genome_length\tSource\tOrganism\tSerotype\tAcronym\tStrain\tNotes\tTitles\n')

for record in SeqIO.parse('initial_download/genbank_files/combined_NCBI_hadv_entries.gb', 'gb'):
    source = record.annotations['source']
    accession = record.name
    bps = len(record.seq)
    title_list = record.annotations['references']
    organism = ''
    serotype = ''
    acronym = ''
    strain = ''
    note = ''

    for feature in record.features:
        feature1 = 'organism'
        feature2 = 'serotype'
        feature3 = 'strain'
        feature4 = 'acronym'
        feature5 = 'note'
        exception1 = 'Human'
        exception2 = 'genotype:'

        try:
            for organism in feature.qualifiers.get(feature1):
                if organism.startswith(exception1):
                    organism = organism
        except Exception:
            pass

        try:
            for serotype in feature.qualifiers.get(feature2):
                serotype = serotype
        except Exception:
            pass

        try:
            for strain in feature.qualifiers.get(feature3):
                strain = strain
        except Exception:
            pass
        try:
            for acronym in feature.qualifiers.get(feature4):
                acronym = acronym
        except Exception:
            pass
        try:
            for note in feature.qualifiers.get(feature5):
                if note.startswith(exception2):
                    note = note
        except Exception:
            pass

    with open(out_file, 'a') as f:
        f.write(
            f'{accession}\t{bps}\t{source}\t{organism}\t{serotype}\t{acronym}\t{strain}\t{note}\t{title_list}\n')
