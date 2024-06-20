import os
import sys
import numpy as np
import pickle
import argparse
import time

def load_distance_matrix(dataset_path):
    if not os.path.exists(dataset_path):
        print(f"File not found: {dataset_path}")
        sys.exit(1)  # Exit the script if the file does not exist
    with open(dataset_path, 'rb') as f:
        return pickle.load(f)

def fast_neighbor_joining(distance_matrix, taxa):
    n = len(taxa)
    clusters = [{i} for i in range(n)]
    cluster_names = taxa.copy()

    while n > 2:
        # Find the two closest clusters
        min_val = np.inf
        for i in range(n):
            for j in range(i + 1, n):
                if distance_matrix[i, j] < min_val:
                    min_val = distance_matrix[i, j]
                    min_pair = (i, j)

        i, j = min_pair

        # Calculate the new distances
        new_distances = [0.5 * (distance_matrix[i, k] + distance_matrix[j, k] - distance_matrix[i, j]) for k in range(n) if k != i and k != j]
        new_row = np.append(new_distances, 0)

        # Update the distance matrix
        mask = np.ones(n, dtype=bool)
        mask[[i, j]] = False
        reduced_matrix = distance_matrix[mask][:, mask]

        new_distance_matrix = np.zeros((n - 1, n - 1))
        new_distance_matrix[:-1, :-1] = reduced_matrix
        new_distance_matrix[-1, :-1] = new_distance_matrix[:-1, -1] = new_row[:-1]

        distance_matrix = new_distance_matrix

        # Update clusters and names
        new_cluster_name = f"({cluster_names[i]}, {cluster_names[j]})"
        clusters.append(clusters[i] | clusters[j])
        cluster_names.append(new_cluster_name)

        for index in sorted([i, j], reverse=True):
            del clusters[index]
            del cluster_names[index]

        n -= 1

    final_distance = distance_matrix[0, 1]
    return f"({cluster_names[0]}:{final_distance / 2},{cluster_names[1]}:{final_distance / 2});"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Fast Neighbor Joining on a dataset.')
    parser.add_argument('--dataset', required=True, help='Path to the dataset file (e.g., data/HIV.npy)')
    parser.add_argument('--output_format', choices=['newick', 'pickle'], default='newick', help='Output format for the tree (default: newick)')
    args = parser.parse_args()

    dist_matrix = load_distance_matrix(args.dataset)
    start_time = time.time()
    taxa = [str(i) for i in range(dist_matrix.shape[0])]
    FNJ_tree = fast_neighbor_joining(dist_matrix, taxa)
    end_time = time.time()
    computation_time = end_time - start_time
    print(f"FNJ computation time: {computation_time} seconds")

    results_dir = 'results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    if args.output_format == 'newick':
        output_path = os.path.join(results_dir, f'FNJ_tree_{os.path.basename(args.dataset)}.newick')
        with open(output_path, 'w') as f:
            f.write(FNJ_tree)
    elif args.output_format == 'pickle':
        output_path = os.path.join(results_dir, f'FNJ_tree_{os.path.basename(args.dataset)}.pickle')
        with open(output_path, 'wb') as f:
            pickle.dump(FNJ_tree, f)

    print(f'FNJ tree saved successfully in {args.output_format} format at {output_path}.')

