import pandas as pd
'''
This script will:
Add a message to missing data fields.
Note: Missing data has been manually confirmed to actually be missing.
    - The specific gene was not annotated
'''

file_path = 'clean_gb_files/logs/HAdV_gene_location_metadata.txt'

df = pd.read_csv(file_path, sep='\t').copy().set_index('Accession_ID')

message1 = 'This genbank entry is missing penton annotations.'
message2 = 'This genbank entry is missing hexon annotations.'
message3 = 'This genbank entry is missing fiber annotations.'

# missing penton locations
df.at['M73260', 'Penton_product'] = message1
df.at['OP270254', 'Penton_product'] = message1
df.at['MG696148', 'Penton_product'] = message1
df.at['MH256657', 'Penton_product'] = message1
df.at['MH256656', 'Penton_product'] = message1
df.at['MH256655', 'Penton_product'] = message1
df.at['MH256654', 'Penton_product'] = message1
df.at['MH256653', 'Penton_product'] = message1

df.at['M73260', 'Penton_location'] = message1
df.at['OP270254', 'Penton_location'] = message1
df.at['MG696148', 'Penton_location'] = message1
df.at['MH256657', 'Penton_location'] = message1
df.at['MH256656', 'Penton_location'] = message1
df.at['MH256655', 'Penton_location'] = message1
df.at['MH256654', 'Penton_location'] = message1
df.at['MH256653', 'Penton_location'] = message1

# missing hexon locations
df.at['M73260', 'Hexon_product'] = message2
df.at['MW816072', 'Hexon_product'] = message2
df.at['MN936178', 'Hexon_product'] = message2
df.at['AF108105', 'Hexon_product'] = message2
df.at['AC_000006', 'Hexon_product'] = message2

df.at['M73260', 'Hexon_location'] = message2
df.at['MW816072', 'Hexon_location'] = message2
df.at['MN936178', 'Hexon_location'] = message2
df.at['AF108105', 'Hexon_location'] = message2
df.at['AC_000006', 'Hexon_location'] = message2

# missing fiber locations
df.at['MN513344', 'Fiber_product'] = message3
df.at['AP014841', 'Fiber_product'] = message3
df.at['M73260', 'Fiber_product'] = message3
df.at['MN936178', 'Fiber_product'] = message3
df.at['KF279629', 'Fiber_product'] = message3
df.at['JN162671', 'Fiber_product'] = message3

df.at['MN513344', 'Fiber_location'] = message3
df.at['AP014841', 'Fiber_location'] = message3
df.at['M73260', 'Fiber_location'] = message3
df.at['MN936178', 'Fiber_location'] = message3
df.at['KF279629', 'Fiber_location'] = message3
df.at['JN162671', 'Fiber_location'] = message3

df.to_csv('clean_gb_files/logs/updated_HAdV_gene_location_metadata.txt', sep='\t')
