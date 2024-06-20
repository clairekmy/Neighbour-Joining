import numpy as np
import os
import subprocess


datasets = [
    'data/E_coli.npy', 'data/Yeast.npy', 'data/HIV.npy', 'data/Aspergillus.npy', 'data/Salmon.npy', 'data/Chicken.npy',
    'data/Mycobacterium_tuberculosis.npy', 'data/Zebrafish.npy', 'data/Corn.npy', 'data/Rice.npy',
    'data/Candida_albicans.npy', 'data/Bacillus_subtilis.npy', 'data/Arabidopsis_thaliana.npy',
    'data/Drosophila_melanogaster.npy', 'data/Human.npy', 'data/Mouse.npy', 'data/Tomato.npy', 'data/Wheat.npy',
    'data/Shrimp.npy', 'data/Streptococcus_pneumoniae.npy'
]

for dataset in datasets:
    data = np.load(dataset)
    print(f"{dataset[:-4]}: {data.shape}")


def generate_pseudo_sequences(num_sequences, sequence_length, output_path):
    # Define the nucleotide alphabet
    nucleotides = ['A', 'C', 'G', 'T']
    
    # Generate random sequences
    pseudo_sequences = [''.join(np.random.choice(nucleotides, sequence_length)) for _ in range(num_sequences)]
    
    # Save the sequences to a FASTA file
    with open(output_path, 'w') as f:
        for i, seq in enumerate(pseudo_sequences):
            f.write(f">seq{i}\n{seq}\n")

# Modify according to your dataset to generate pseudo sequences for various organisms
generate_pseudo_sequences(350, 2936, 'data/E_coli_pseudo.fasta')
generate_pseudo_sequences(350, 2467, 'data/Yeast_pseudo.fasta')
generate_pseudo_sequences(350, 4835, 'data/Arabidopsis_thaliana_pseudo.fasta')
generate_pseudo_sequences(350, 999, 'data/Drosophila_melanogaster_pseudo.fasta')
generate_pseudo_sequences(593, 4884, 'data/Salmon_pseudo.fasta')
generate_pseudo_sequences(1200, 1498, 'data/Mouse_pseudo.fasta')
generate_pseudo_sequences(1494, 1000, 'data/Tomato_pseudo.fasta')
generate_pseudo_sequences(1200, 5000, 'data/HIV_pseudo.fasta')
generate_pseudo_sequences(1500, 3000, 'data/Zebrafish_pseudo.fasta')
generate_pseudo_sequences(1199, 1000, 'data/Wheat_pseudo.fasta')
generate_pseudo_sequences(1493, 4953, 'data/Aspergillus_pseudo.fasta')
generate_pseudo_sequences(1200, 4476, 'data/Streptococcus_pneumoniae_pseudo.fasta')
generate_pseudo_sequences(1500, 4859, 'data/Rice_pseudo.fasta')
generate_pseudo_sequences(600, 4486, 'data/Chicken_pseudo.fasta')
generate_pseudo_sequences(600, 4996, 'data/Mycobacterium_tuberculosis_pseudo.fasta')
generate_pseudo_sequences(600, 4981, 'data/Shrimp_pseudo.fasta')
generate_pseudo_sequences(880, 4269, 'data/Candida_albicans_pseudo.fasta')
generate_pseudo_sequences(900, 4444, 'data/Bacillus_subtilis_pseudo.fasta')
generate_pseudo_sequences(893, 4694, 'data/Corn_pseudo.fasta')
generate_pseudo_sequences(900, 1500, 'data/Human_pseudo.fasta')







# Add the FastTree directory to the PATH
os.environ['PATH'] += os.pathsep + '/Users/xxx/opt/anaconda3/bin'




def run_fasttree(fasta_file, output_tree):
    fasttree_cmd = f"fasttree -nt < {fasta_file} > {output_tree}"
    subprocess.run(fasttree_cmd, shell=True, check=True)

# Run FastTree on the pseudo sequences for selected organisms
run_fasttree('data/E_coli_pseudo.fasta', 'Reference/E_coli.tree'),
run_fasttree('data/Yeast_pseudo.fasta', 'Reference/Yeast.tree'),
run_fasttree('data/HIV_pseudo.fasta', 'Reference/HIV.tree'),
run_fasttree('data/Aspergillus_pseudo.fasta', 'Reference/Aspergillus.tree'),
run_fasttree('data/Salmon_pseudo.fasta', 'Reference/Salmon.tree'),
run_fasttree('data/Chicken_pseudo.fasta', 'Reference/Chicken.tree'),
run_fasttree('data/Mycobacterium_tuberculosis_pseudo.fasta', 'Reference/Mycobacterium_tuberculosis.tree'),
run_fasttree('data/Zebrafish_pseudo.fasta', 'Reference/Zebrafish.tree'),
run_fasttree('data/Corn_pseudo.fasta', 'Reference/Corn.tree'),
run_fasttree('data/Rice_pseudo.fasta', 'Reference/Rice.tree'),
run_fasttree('data/Candida_albicans_pseudo.fasta', 'Reference/Candida_albicans.tree'),
run_fasttree('data/Bacillus_subtilis_pseudo.fasta', 'Reference/Bacillus_subtilis.tree'),
run_fasttree('data/Arabidopsis_thaliana_pseudo.fasta', 'Reference/Arabidopsis_thaliana.tree'),
run_fasttree('data/Drosophila_melanogaster_pseudo.fasta', 'Reference/Drosophila_melanogaster.tree'),
run_fasttree('data/Human_pseudo.fasta', 'Reference/Human.tree'),
run_fasttree('data/Mouse_pseudo.fasta', 'Reference/Mouse.tree'),
run_fasttree('data/Tomato_pseudo.fasta', 'Reference/Tomato.tree'),
run_fasttree('data/Wheat_pseudo.fasta', 'Reference/Wheat.tree'),
run_fasttree('data/Shrimp_pseudo.fasta', 'Reference/Shrimp.tree'),
run_fasttree('data/Streptococcus_pneumoniae_pseudo.fasta', 'Reference/Streptococcus_pneumoniae.tree')

