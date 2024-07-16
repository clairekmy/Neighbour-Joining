import os
import sys
import pickle
import subprocess
import time
import numpy as np
import argparse
import matplotlib.pyplot as plt
import utils  

def save_distance_matrix(dataset_path):
    # Load the dataset
    data = np.load(dataset_path)

    n_leaves = data.shape[0]
    seq_length = data.shape[1]

    # Convert the sequences to one-hot encoding
    data_onehot = utils.convert_data_str_to_onehot(data)

    # Calculate the distance matrix
    way_1_tic = time.process_time()
    dist_matrix = np.zeros([n_leaves, n_leaves])
    for ki in range(n_leaves):
        dist_matrix[ki, ki+1:] += np.abs(data_onehot[ki+1:, :, :] - data_onehot[ki, :, :]).sum(axis=-1).sum(axis=-1)
    dist_matrix += dist_matrix.T
    dist_matrix = dist_matrix / (2 * data.shape[1])

    # Adjust the distance matrix for the Jukes-Cantor model
    epsilon = 1e-10
    dist_matrix = np.clip(dist_matrix, 0, 1 - 4 * epsilon / 3)
    dist_matrix = (-3 / 4) * np.log(1 - dist_matrix * 4 / 3 + epsilon)
    dist_matrix = (dist_matrix + dist_matrix.T) / 2
    dist_matrix = np.nan_to_num(dist_matrix, nan=0.0)
    np.fill_diagonal(dist_matrix, 0)
    way_1_toc = time.process_time()

    results_dir = 'results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Save the distance matrix with a dataset-specific filename
    dist_matrix_output_path = os.path.join(results_dir, f'distance_matrix_{os.path.basename(dataset_path)}.pickle')
    with open(dist_matrix_output_path, 'wb') as f:
        pickle.dump(dist_matrix, f)

    # Print the timings
    print('*******************')
    print('DONE=>>>>> ' + dataset_path)
    print('Distance matrix calculation time=', way_1_toc - way_1_tic)
    print('*******************')
    
    return dist_matrix, way_1_toc - way_1_tic

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run the analysis on given datasets.')
    parser.add_argument('--datasets', nargs='+', default=[
        'Enterobacter_cloacae', 'Drosophila_melanogaster', 'Shigella_flexneri', 'Arabidopsis_thaliana', 'Helicobacter_pylori',
        'E_coli', 'Pseudomonas_aeruginosa', 'Yersinia_pestis', 'Yeast', 'Legionella_pneumophila',
        'Salmon', 'Vibrio_cholerae', 'Bartonella_henselae', "Baker's_yeast", 'Chicken', 'Shrimp',
        'Treponema_pallidum', 'Listeria_monocytogenes', 'Mycobacterium_tuberculosis', 'Francisella_tularensis',
        'Candida_albicans', 'Human', 'Corn', 'Staphylococcus_aureus', 'Bacillus_subtilis', 'Vibrio_vulnificus',
        'Lactobacillus_casei', 'Caenorhabditis_elegans', 'Corynebacterium_diphtheriae', 'Streptococcus_pneumoniae',
        'Wheat', 'HIV', 'Mouse'],
        help='Names of the datasets')
    parser.add_argument('--num_of_seed', type=int, default=10, help='Number of seeds (e.g., 5, 10, 20)')

    args = parser.parse_known_args()[0]

    running_times = {'NJ': [], 'SparseNJ': [], 'FNJ': [], 'RNJ': []}

    for dataset in args.datasets:
        sparse_dataset_path = os.path.join('data', f'{dataset}.npy')
        results_dataset_path = os.path.join('results', f'distance_matrix_{dataset}.npy.pickle')

        # Calculate the distance matrix and get the calculation time
        dist_matrix, dist_matrix_time = save_distance_matrix(sparse_dataset_path)

        # Save the distance matrix for reuse
        with open(results_dataset_path, 'wb') as f:
            pickle.dump(dist_matrix, f)

        # Run NJ first and include the distance matrix calculation time
        tic = time.perf_counter()
        try:
            result = subprocess.run([sys.executable, 'NJ.py', '--dataset', dataset],
                                    capture_output=True, text=True, check=True)
            toc = time.perf_counter()
            elapsed_time = toc - tic + dist_matrix_time  # Add the distance matrix time

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
                elapsed_time = toc - tic + (dist_matrix_time if algorithm in ['FNJ', 'RNJ'] else 0)

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
