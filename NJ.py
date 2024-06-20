import os
import sys
import pickle
import time
import numpy as np
import argparse
from skbio import DistanceMatrix
from skbio.tree import nj

import sys

def load_distance_matrix(dataset_path):
    if not os.path.exists(dataset_path):
        print(f"File not found: {dataset_path}")
        sys.exit(1)  # Exit the script if the file does not exist
    with open(dataset_path, 'rb') as f:
        dist_matrix = pickle.load(f)
    return dist_matrix


def process_dataset(dataset_name):
    """Process a single dataset to generate and save a Neighbor-Joining tree."""
    dist_matrix_path = os.path.join('results', f'distance_matrix_{dataset_name}.npy.pickle')
    if not os.path.exists(dist_matrix_path):
        print(f"Distance matrix file not found: {dist_matrix_path}")
        return

    dist_matrix = load_distance_matrix(dist_matrix_path)
    taxa = [str(i) for i in range(len(dist_matrix))]

    # Calculate NJ tree
    NJ_tic = time.perf_counter()
    dm = DistanceMatrix(dist_matrix, taxa)
    nj_tree = nj(dm)
    NJ_toc = time.perf_counter()

    # Save the NJ tree using pickle
    output_path = os.path.join('results', f'NJ_tree_{dataset_name}.pickle')
    try:
        with open(output_path, 'wb') as f:
            pickle.dump(nj_tree, f)
        print(f"NJ tree saved successfully to {output_path}")
    except Exception as e:
        print(f"Error saving NJ tree: {e}")

    print('*******************')
    print(f'DONE=>>>>> {dataset_name}')
    print(f'NJ computation time: {NJ_toc - NJ_tic} seconds')
    print('*******************')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', required=True, help='Name of the dataset (e.g., HIV)')
    args = parser.parse_args()

    process_dataset(args.dataset)
