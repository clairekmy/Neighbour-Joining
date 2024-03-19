import numpy as np
import time
import pickle
import utils  
import heapq
from tqdm import tqdm
import os
import argparse


def fast_neighbor_joining(distance_matrix, taxa):
    n = len(distance_matrix)
    clusters = [{i} for i in range(n)]
    cluster_names = taxa.copy()
    active_clusters = set(range(n))

    # Initialize a priority queue with distances and cluster pairs
    pq = []
    for i in range(n):
        for j in range(i + 1, n):
            heapq.heappush(pq, (distance_matrix[i, j], (i, j)))

    with tqdm(total=n - 2, desc="Merging clusters") as pbar:
        while len(active_clusters) > 2:
            # Find the pair of active clusters with the minimum distance
            while True:
                dist, (i, j) = heapq.heappop(pq)
                if i in active_clusters and j in active_clusters:
                    break

            # Merge the clusters
            new_cluster = clusters[i] | clusters[j]
            new_cluster_name = f"({cluster_names[i]}:{distance_matrix[i, j]/2},{cluster_names[j]}:{distance_matrix[i, j]/2})"
            clusters.append(new_cluster)
            cluster_names.append(new_cluster_name)
            k = len(clusters) - 1  # This is the new cluster's index
            active_clusters.add(k)

            # Update distances to the new cluster
            new_distances = np.zeros(n + 1)
            for l in active_clusters:
                if l != k:
                    d = (distance_matrix[i, l] + distance_matrix[j, l] - distance_matrix[i, j]) / 2
                    new_distances[l] = d
                    heapq.heappush(pq, (d, (k, l)))

            # Update the distance matrix for the new cluster
            new_matrix = np.zeros((n + 1, n + 1))
            new_matrix[:-1, :-1] = distance_matrix  # Copy the old distances
            new_matrix[-1, :-1] = new_distances[:-1]  # Update with new distances
            new_matrix[:-1, -1] = new_distances[:-1]  # Symmetric update
            distance_matrix = new_matrix

            # Remove the merged clusters from active clusters and adjust for new indexing
            active_clusters.remove(i)
            active_clusters.remove(j)

            pbar.update(1)  # Update the progress bar

            # Adjust the size of n to reflect the new distance matrix size
            n = len(distance_matrix)

    # Merge the last two clusters
    i, j = active_clusters
    new_cluster_name = f"({cluster_names[i]}:{distance_matrix[i, j]/2},{cluster_names[j]}:{distance_matrix[i, j]/2});"
    return new_cluster_name





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Fast Neighbor Joining on a dataset.')
    parser.add_argument('--dataset', required=True, help='Path to the dataset file (e.g., data/HIV.npy)')
    args = parser.parse_args()

    # Load the dataset
    data = np.load(args.dataset)

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

    # Construct the FNJ tree
    FNJ_tic = time.process_time()
    taxa = [str(i) for i in range(n_leaves)]
    FNJ_tree = fast_neighbor_joining(dist_matrix, taxa)
    FNJ_toc = time.process_time()

    results_dir = 'results'
    

    # Save the FNJ tree with a dataset-specific filename
    fnj_output_path = os.path.join(results_dir, f'FNJ_tree_{args.dataset}.pickle')
    with open(fnj_output_path, 'wb') as f:
        pickle.dump(FNJ_tree, f)


    # Print the timings
    print('*******************')
    print('DONE=>>>>> ' + args.dataset)
    print('FNJ clustering time=', FNJ_toc - FNJ_tic)
    print('*******************')
