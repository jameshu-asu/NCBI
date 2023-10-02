from Bio import Entrez

'''By James C. Hu
This script will:
1) Search for HAdV genbank files from specified NCBI databases using search terms.
2) Combine all like files into a single genbank file.
    * Search terms is a single comma-seperated string of ids.
3) Generate a log file of how many entries were returned for each HAdV type.
'''

#____NCBI_identification____#
Entrez.email = 'email@domain.com'
handle = Entrez.einfo()
NCBI_db_list = Entrez.read(handle)
print('Below is a list of NCBI databases:')
print(NCBI_db_list)
print('\n')

#_____Variables/File_names/Paths____#
# count start/end double as HAdV types. As of 8.13.23 there are 113 differnt HAdV genotypes
COUNT_START = 1
COUNT_END = 114
log_file = 'initial_download/download_logs/NCBI_total_entries_by_type_log.txt'
gb_file = 'initial_download/genbank_files/combined_NCBI_hadv_entries.gb'


#____Create_directory_structure____#
os.makedirs('initial_download/genbank_files')
os.makedirs('initial_download/download_logs')

#____Instantiate_out_files____#
# The file below will display the specific HAdV type along with the total record count returned for that type.
with open(log_file, 'w') as f:
    f.write('HAdV_type\tTotal_recores\n')

# The file below will display all entries retrieved from NCBI
with open(gb_file, 'w') as f:
    pass

#____Accessing_NCBI_DB____#
# Using count and a while loop to automate search/download of all 113 HAdV types
while COUNT_START < COUNT_END:
    print(f'HAdV type {COUNT_START}')
    handle = Entrez.esearch(db='nuccore',
                            term=f'human adenovirus {COUNT_START} complete genome')
    records = Entrez.read(handle)
    total_records = records['Count']
    print(f'Total records: {total_records}\n')
    # append record data to log file
    with open(log_file, 'a') as f:
        f.write(f'{COUNT_START}\t{total_records}\n')
    handle = Entrez.esearch(db='nuccore',
                            term=f'human adenovirus {COUNT_START} complete genome',
                            retmax=total_records)
    records = Entrez.read(handle)
    # Generate single search term: id_str
    id_list = records['IdList']
    id_str = ','.join(id_list)
    print('Combined NCBI nuccore internal id list:')
    print(f'{id_str}\n')
    handle = Entrez.efetch(db='nuccore',
                           id=id_str,
                           rettype='gb',
                           retmode='text')
    text = handle.read()
    # Search restuls will be written to a single gb file.
    with open(gb_file, 'a') as f:
        f.write(text)
    COUNT_START += 1
