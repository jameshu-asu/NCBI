from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_consensus
from Bio.Phylo.Consensus import majority_consensus
import matplotlib.pyplot as plt
import os

'''By James C. Hu
This script will:
1) Generate png and newick tree file from MSA aligned fasta file.
2) Generate bootstrapped png and newick tree file from MSA aligned fasta file.
'''


def tree_builder_fa(alined_fasta_file: str, nwk_outfile: str) -> None:
    '''Builds Phylogenetic tree using alignment file generated during MSA
    input: Aligned fasta file
    output: Phylogenetic tree png, Newick format tree file.

    alined_fasta_file: fasta infile
    nwk_outfile: name of nwk outfile
    '''
    aligned_fasta = AlignIO.read(alined_fasta_file, 'fasta')
    print(aligned_fasta)
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(aligned_fasta)
    print('\nDistance Matrix\n==========================')
    print(distance_matrix)
    constructor = DistanceTreeConstructor()
    # options: nj (neighbor-joining) or upgma (Unweighted Pair Group Method with Arithemtic Mean)
    tree = constructor.nj(distance_matrix)
    print('\nPhylogenetic Tree\n==========================')
    Phylo.draw_ascii(tree)
    Phylo.write(tree, nwk_outfile, 'newick')
    plt.rc('font', size=5)
    fig, ax = plt.subplots(figsize=(7, 9))
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.savefig('tree.png', dpi=1080)
    # plt.show()
    return None


def bootstrap_consensus_tree(aligned_fasta_file: str, replicates: int) -> None:
    '''Builds consensus tree based on n bootstrapped replicate trees.
    input: aligned fasta
    output: bootstrapped n replicate consensus tree.
    '''
    aligned_fasta = AlignIO.read(aligned_fasta_file, 'fasta')
    calculator = DistanceCalculator('blosum62')
    constructor = DistanceTreeConstructor(calculator)
    tree = bootstrap_consensus(aligned_fasta, replicates, constructor, majority_consensus)
    Phylo.write(tree, 'consensus_tree.nwk', 'newick')
    plt.rc('font', size=5)
    fig, ax = plt.subplots(figsize=(7, 9))
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.savefig('bootstrap_consensus_tree.png', dpi=1080)
    Phylo.draw_ascii(tree)
    return None


os.makedirs('clean_gb_files/combined_fastas/trees')

in_path = 'clean_gb_files/combined_fastas/MSA_output'
out_path = 'clean_gb_files/combined_fastas/trees'

aligned_penton_fa = os.path.join(in_path, 'MSA_aligned_penton_output.fasta')
aligned_hexon_fa = os.path.join(in_path, 'MSA_aligned_hexon_output.fasta')
aligned_fiber_fa = os.path.join(in_path, 'MSA_aligned_fiber_output.fasta')

aligned_penton_tree = os.path.join(out_path, 'aligned_penton_tree.nwk')
aligned_hexon_tree = os.path.join(out_path, 'aligned_hexon_tree.nwk')
aligned_fiber_tree = os.path.join(out_path, 'aligned_fiber_tree.nwk')

tree_builder_fa(aligned_penton_fa, aligned_penton_tree)
tree_builder_fa(aligned_hexon_fa, aligned_hexon_tree)
tree_builder_fa(aligned_fiber_fa, aligned_fiber_tree)
