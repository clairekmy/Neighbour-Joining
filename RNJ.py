import numpy as np
import time
import pickle
import os
import argparse

def load_distance_matrix(dataset_path):
    if not os.path.exists(dataset_path):
        print(f"File not found: {dataset_path}")
        sys.exit(1)  # Exit the script if the file does not exist
    with open(dataset_path, 'rb') as f:
        return pickle.load(f)

def relaxed_neighbor_joining(distance_matrix, taxa):
    n = len(taxa)
    clusters = [{i} for i in range(n)]
    cluster_names = taxa.copy()

    while n > 2:
        min_val = np.inf
        for x in range(n):
            for y in range(x + 1, n):
                if distance_matrix[x, y] < min_val:
                    min_val = distance_matrix[x, y]
                    i, j = x, y

        new_distances = [(distance_matrix[i, k] + distance_matrix[j, k] - distance_matrix[i, j]) / 2 for k in range(n) if k != i and k != j]
        new_distances.append(0)  # Distance to itself

        new_distance_matrix = np.zeros((n - 1, n - 1))
        for x in range(n - 2):
            for y in range(n - 2):
                xi, yi = (x if x < i else x + 1), (y if y < j else y + 1)
                new_distance_matrix[x, y] = distance_matrix[xi, yi]
        for x in range(n - 1):
            new_distance_matrix[x, -1] = new_distance_matrix[-1, x] = new_distances[x]

        distance_matrix = new_distance_matrix
        new_cluster_name = f"({cluster_names[i]}:{min_val / 2},{cluster_names[j]}:{min_val / 2})"
        clusters.append(clusters[i] | clusters[j])
        cluster_names.append(new_cluster_name)
        for index in sorted([i, j], reverse=True):
            del clusters[index]
            del cluster_names[index]

        n -= 1

    final_distance = distance_matrix[0, 1]
    return f"({cluster_names[0]}:{final_distance / 2},{cluster_names[1]}:{final_distance / 2});"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Relaxed Neighbor Joining on a dataset.')
    parser.add_argument('--dataset', required=True, help='Path to the dataset file (e.g., data/HIV.npy)')
    args = parser.parse_args()

    dist_matrix = load_distance_matrix(args.dataset)
    start_time = time.time()
    taxa = [str(i) for i in range(dist_matrix.shape[0])]
    RNJ_tree = relaxed_neighbor_joining(dist_matrix, taxa)
    end_time = time.time()
    computation_time = end_time - start_time
    print(f"RNJ computation time: {computation_time} seconds")

    results_dir = 'results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    rnj_output_path = os.path.join(results_dir, f'RNJ_tree_{os.path.basename(args.dataset)}.newick')
    with open(rnj_output_path, 'w') as f:
        f.write(RNJ_tree)

    # Optionally save computation metadata in a pickle file
    metadata = {
        'computation_time': computation_time,
        'num_taxa': len(taxa),
        'algorithm': 'Relaxed Neighbor Joining'
    }
    metadata_path = os.path.join(results_dir, f'RNJ_metadata_{os.path.basename(args.dataset)}.pickle')
    with open(metadata_path, 'wb') as f:
        pickle.dump(metadata, f)

    print(f'RNJ tree saved successfully to {rnj_output_path}')
    print(f'Metadata saved successfully to {metadata_path}')
