import os
import sys
import subprocess
import time
import numpy as np
import argparse
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def load_npy(file_path):
    return np.load(file_path)

def npy_to_fasta(npy_file, fasta_file):
    sequences = np.load(npy_file)
    seq_list = [''.join(seq) for seq in sequences.astype(str)]
    records = []
    for i, seq in enumerate(seq_list):
        record = SeqRecord(Seq(seq), id=f"seq_{i+1}", description="")
        records.append(record)
    with open(fasta_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')
    print(f"Converted {npy_file} to {fasta_file}")

def construct_tree_with_clearcut(fasta_file):
    result = subprocess.run(['/Users/kangmingyue/NUS/Research/Phylogenetic_Massive/clearcut/clearcut', '--in', fasta_file, '--out', 'temp_tree.newick', '--alignment', '--DNA'], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"clearcut failed: {result.stderr}")
    with open('temp_tree.newick', 'r') as f:
        newick_tree = f.read().strip()
    os.remove('temp_tree.newick')
    return newick_tree

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run the analysis on given datasets.')
    parser.add_argument('--datasets', nargs='+', default=[
        'E_coli', 'Yeast', 'Arabidopsis_thaliana', 'Drosophila_melanogaster', 'Salmon',
        'Chicken', 'Mycobacterium_tuberculosis', 'Shrimp', 'Candida_albicans', 'Corn',
        'Bacillus_subtilis', 'Human', 'Mouse', 'HIV', 'Streptococcus_pneumoniae',
        'Wheat', 'Aspergillus', 'Tomato', 'Zebrafish', 'Rice'],
        help='Names of the datasets')
    parser.add_argument('--num_of_seed', type=int, default=10, help='Number of seeds (e.g., 5, 10, 20)')

    args = parser.parse_known_args()[0]

    running_times = {'NJ': [], 'SparseNJ': [], 'FNJ': [], 'RNJ': []}

    for dataset in args.datasets:
        npy_dataset_path = os.path.join('data', f'{dataset}.npy')
        fasta_dataset_path = os.path.join('data', f'{dataset}_aligned.fasta')

        if not os.path.exists(npy_dataset_path):
            print(f"File not found: {npy_dataset_path}")
            continue

        # Convert NPY to FASTA for RNJ or clearcut
        npy_to_fasta(npy_dataset_path, fasta_dataset_path)

        # Load the NPY file
        data = load_npy(npy_dataset_path)

        # Run NJ first
        tic = time.perf_counter()
        try:
            result = subprocess.run([sys.executable, 'NJ.py', '--dataset', npy_dataset_path],
                                    capture_output=True, text=True, check=True)
            toc = time.perf_counter()
            elapsed_time = toc - tic

            running_times['NJ'].append(elapsed_time)

            print('*******************')
            print(f'DONE=>>>>> {dataset} with NJ')
            print(f'NJ computation time: {elapsed_time} seconds')
            print('*******************')

        except subprocess.CalledProcessError as e:
            print(f"Error running NJ.py: {e}")
            print(f"Stderr: {e.stderr}")
            continue

        # Run other algorithms with NPY files
        for algorithm in ['SparseNJ', 'FNJ']:
            tic = time.perf_counter()
            try:
                result = subprocess.run([sys.executable, f'{algorithm}.py', '--dataset', npy_dataset_path],
                                        capture_output=True, text=True, check=True)
                toc = time.perf_counter()
                elapsed_time = toc - tic

                running_times[algorithm].append(elapsed_time)

                print('*******************')
                print(f'DONE=>>>>> {dataset} with {algorithm}')
                print(f'{algorithm} computation time: {elapsed_time} seconds')
                print('*******************')

            except subprocess.CalledProcessError as e:
                print(f"Error running {algorithm}.py: {e}")
                print(f"Stderr: {e.stderr}")

        # Run RNJ with clearcut
        tic = time.perf_counter()
        try:
            RNJ_tree = construct_tree_with_clearcut(fasta_dataset_path)
            toc = time.perf_counter()
            elapsed_time = toc - tic

            running_times['RNJ'].append(elapsed_time)

            print('*******************')
            print(f'DONE=>>>>> {dataset} with RNJ')
            print(f'RNJ computation time: {elapsed_time} seconds')
            print('*******************')

        except RuntimeError as e:
            print(f"Error running RNJ with clearcut: {e}")

    plt.figure(figsize=(12, 6))
    markers = ['o', 's', 'D', 'x']
    colors = ['#1f77b4', '#2ca02c', '#d62728', '#ff7f0e']
    labels = ['NJ', 'SparseNJ', 'FNJ', 'RNJ']

    print("Running times:", running_times)

    for label, marker, color in zip(labels, markers, colors):
        if running_times[label]:
            plt.scatter(args.datasets, running_times[label], s=20, marker=marker, label=label, color=color, edgecolor='black')
            plt.plot(args.datasets, running_times[label], color=color, linestyle='-', linewidth=2)
        else:
            print(f"No running times for {label}")

    plt.xlabel('Dataset', fontsize=16, labelpad=12)
    plt.ylabel('Elapsed Time (seconds)', fontsize=16, labelpad=12)
    plt.title('Elapsed Time Comparison', fontsize=18, pad=22)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14, frameon=True, shadow=True)
    plt.tight_layout()
    plt.savefig('elapsed_time_comparison_scatter_all.png', dpi=300, bbox_inches='tight')
    plt.show()
