from Bio.Align.Applications import MafftCommandline
import os

'''By James C. Hu
This script will:
1) Make a directory to house MSA outputs
2) Generate MSA (MAFFT) output for combined files.
'''


def MSA_mafft(infile: str, output: str) -> None:
    '''This funciton perfrom MSA on combined consensus files using MAFFT
    input: Combined consenssus fasta file.
    output: aligned.fasta
    '''
    mafft_exe = '/Users/account/anaconda3/bin/mafft'
    mafft_cline = MafftCommandline(mafft_exe, input=infile)
    stdout, stderr = mafft_cline()
    with open(f'{output}', 'w') as handle:
        handle.write(stdout.upper())
    return None


os.makedirs('clean_gb_files/combined_fastas/MSA_output')

in_path = 'clean_gb_files/combined_fastas'
combined_penton_fa_path = os.path.join(in_path, 'combined_penton_sequences.fasta')
combined_hexon_fa_path = os.path.join(in_path, 'combined_hexon_sequences.fasta')
combined_fiber_fa_path = os.path.join(in_path, 'combined_fiber_sequences.fasta')

out_path = 'clean_gb_files/combined_fastas/MSA_output'
aligned_penton_fa = os.path.join(out_path, 'MSA_aligned_penton_output.fasta')
aligned_hexon_fa = os.path.join(out_path, 'MSA_aligned_hexon_output.fasta')
aligned_fiber_fa = os.path.join(out_path, 'MSA_aligned_fiber_output.fasta')


MSA_mafft(combined_penton_fa_path, aligned_penton_fa)
MSA_mafft(combined_hexon_fa_path, aligned_hexon_fa)
MSA_mafft(combined_fiber_fa_path, aligned_fiber_fa)
