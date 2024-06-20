import os
import sys
import pickle
import subprocess
import time
import numpy as np
import argparse
import matplotlib.pyplot as plt

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
        sparse_dataset_path = os.path.join('data', f'{dataset}.npy')
        results_dataset_path = os.path.join('results', f'distance_matrix_{dataset}.npy.pickle')

        if os.path.exists(results_dataset_path):
            with open(results_dataset_path, 'rb') as f:
                data_results = pickle.load(f)
            print(f"Distance matrix for {dataset} loaded for NJ, RNJ, FNJ.")
        else:
            print(f"File not found: {results_dataset_path}")
            continue

        # Run NJ first
        tic = time.perf_counter()
        try:
            result = subprocess.run([sys.executable, 'NJ.py', '--dataset', dataset],
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

        for algorithm in ['SparseNJ', 'FNJ', 'RNJ']:
            tic = time.perf_counter()
            try:
                dataset_path = sparse_dataset_path if algorithm == 'SparseNJ' else results_dataset_path

                result = subprocess.run([sys.executable, f'{algorithm}.py', '--dataset', dataset_path],
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

    plt.figure(figsize=(12, 6))
    markers = ['o', 's', 'D', 'x']
    colors = ['#1f77b4', '#2ca02c', '#d62728', '#ff7f0e']
    labels = ['NJ', 'SparseNJ', 'FNJ', 'RNJ']

    print("Running times:", running_times)

    for label, marker, color in zip(labels, markers, colors):
        if running_times[label]:
            plt.scatter(args.datasets, running_times[label], s=20, marker=marker, label=label, color=color, edgecolor='black')
            plt.plot(args.datasets, running_times[label], color=color, linestyle='-', linewidth=2)  # Line connecting the points
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
