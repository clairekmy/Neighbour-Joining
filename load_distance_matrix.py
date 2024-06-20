import numpy as np
import os
import pickle
import utils
import time
import argparse

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

    # Save the distance matrix with a dataset-specific filename
    dist_matrix_output_path = os.path.join(results_dir, f'distance_matrix_{os.path.basename(dataset_path)}.pickle')
    with open(dist_matrix_output_path, 'wb') as f:
        pickle.dump(dist_matrix, f)

    # Print the timings
    print('*******************')
    print('DONE=>>>>> ' + dataset_path)
    print('Distance matrix calculation time=', way_1_toc - way_1_tic)
    print('*******************')



datasets = [
    'E_coli', 'Yeast', 'Arabidopsis_thaliana', 'Drosophila_melanogaster', 'Salmon',
    'Chicken', 'Mycobacterium_tuberculosis', 'Shrimp', 'Candida_albicans', 'Corn',
    'Bacillus_subtilis', 'Human', 'Mouse', 'HIV', 'Streptococcus_pneumoniae',
    'Wheat', 'Aspergillus', 'Tomato', 'Zebrafish', 'Rice'
]




# Process each dataset
for dataset in datasets:
    # Construct the full path to the dataset file
    dataset_path = os.path.join('data', f'{dataset}.npy')
    
    # Call the save_distance_matrix function for each dataset
    save_distance_matrix(dataset_path)



