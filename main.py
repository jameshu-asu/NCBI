import subprocess
import time
import sys

'''
Workflow script for accessing, qcing and parsing HAdV gene sequences from combined NCBI gb file.
'''


def time_subprocess(timer_log: str, script: str) -> None:
    with open(timer_log, 'a') as f:
        start = time.time()
        subprocess.call(['python', script])
        end = time.time()
        delta = end - start
        f.write(f'{script}\t{delta}\n')
    return None


timer_log = 'timer_log.txt'

with open(timer_log, 'w') as f:
    f.write(f'Subprocess\tElapsed_time\n')

# script 1: downloading gb files from NCBI db
print(f'Begining downloading HAdV genbank entries from NCBI nucleotide database.\n')
print(f'Running downloading_files_1.py\n')
try:
    time_subprocess(timer_log, 'download_files_1.py')
    print(f'Download complete.\n')
except Exception as e:
    print('Exiting workflow due to the following error:')
    print(e)
    sys.exit()


# script 2: generate metadata file for entries downloaded using script 1
print(f'Begining generating download metadata file.')
print(f'Running generate_gb_metadata_file_2.py\n')
try:
    time_subprocess(timer_log, 'generate_gb_metadata_file_2.py')
    print(f'Download metadata file generated.\n')
except Exception as e:
    print('Exiting workflow due to the following error:')
    print(e)
    sys.exit()


# script 3: qc metadata file generated from script 2
print(f'Begining qc step1 on metadata file.')
print(f'Running qc_gb_metadata_file_3.py\n')
try:
    time_subprocess(timer_log, 'qc_gb_metadata_file_3.py')
    print(f'QC step1 complete.\n')
except Exception as e:
    print('Exiting workflow due to the following error:')
    print(e)
    sys.exit()

# script 4: sorting missmatched samples
print(f'Begining qc step 2 on metadata file.')
print(f'Running missmatch_sort_4.py\n')
try:
    time_subprocess(timer_log, 'missmatch_sort_4.py')
    print(f'QC step2 complete.\n')
except Exception as e:
    print('Exiting workflow due to the following error:')
    print(e)
    sys.exit()

# script 5: generate
print(f'Parsing gb file.')
print(f'Running combine_hadv_entries_by_type_5.py\n')
try:
    time_subprocess(timer_log, 'combine_hadv_entries_by_type_5.py')
    print(f'Accessions have been split by type.\n')
except Exception as e:
    print('Exiting workflow due to the following error:')
    print(e)
    sys.exit()

# script 6: generate
print(f'Generating feature list.')
print(f'Running pull_entry_features_6.py\n')
try:
    time_subprocess(timer_log, 'pull_entry_features_6.py')
    print(f'Feature list created.\n')
except Exception as e:
    print('Exiting workflow due to the following error:')
    print(e)
    sys.exit()

# script 7: generate
print(f'Generating gene location list.')
print(f'Running clean_gene_location_metadata_7.py\n')
try:
    time_subprocess(timer_log, 'clean_gene_location_metadata_7.py')
    print(f'Gene location metadata file created.\n')
except Exception as e:
    print('Exiting workflow due to the following error:')
    print(e)
    sys.exit()

# script 8: generate
print(f'Moving gene fasta files.')
print(f'Running sort_gene_gb_files_8.py\n')
try:
    time_subprocess(timer_log, 'sort_gene_gb_files_8.py')
    print(f'All processes complete.\n')
except Exception as e:
    print('Exiting workflow due to the following error:')
    print(e)
    sys.exit()
