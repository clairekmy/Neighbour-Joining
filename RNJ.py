import numpy as np
import time
import pickle
import utils  
from tqdm import tqdm

from tqdm import tqdm


def relaxed_neighbor_joining(distance_matrix, taxa):
    clusters = [{i} for i in range(len(taxa))]
    cluster_names = taxa.copy()
    total_iterations = len(taxa) - 2  # Total number of iterations required

    with tqdm(total=total_iterations, desc="Merging clusters") as pbar:
        while len(clusters) > 2:
            n = len(distance_matrix)
            total_distances = np.sum(distance_matrix, axis=1)
            q_matrix = (n - 2) * distance_matrix - total_distances[:, None] - total_distances[None, :]

            i, j = np.unravel_index(np.argmin(q_matrix), q_matrix.shape)

            new_cluster_distances = (distance_matrix[i, :] + distance_matrix[j, :] - distance_matrix[i, j]) / 2

            new_cluster = clusters[i] | clusters[j]
            new_cluster_name = f"({cluster_names[i]}, {cluster_names[j]})"

            clusters.append(new_cluster)
            cluster_names.append(new_cluster_name)

            # Create a new distance matrix
            new_distance_matrix = np.zeros((n - 1, n - 1))
            for a in range(n - 1):
                for b in range(a + 1, n - 1):
                    ai = a if a < i else a + 1
                    bi = b if b < j else b + 1
                    new_distance_matrix[a, b] = new_distance_matrix[b, a] = distance_matrix[ai, bi]
            new_distance_matrix[-1, :-1] = new_cluster_distances[:-2]  # Use elements up to the second last
            new_distance_matrix[:-1, -1] = new_cluster_distances[:-2]  # Use elements up to the second last

            distance_matrix = new_distance_matrix

            if i < j:
                del clusters[j]
                del clusters[i]
                del cluster_names[j]
                del cluster_names[i]
            else:
                del clusters[i]
                del clusters[j]
                del cluster_names[i]
                del cluster_names[j]

            pbar.update(1)  # Update the progress bar

    i, j = 0, 1
    new_cluster_name = f"({cluster_names[i]}, {cluster_names[j]})"
    return new_cluster_name




# Load the sequences
data = np.load('hiv_sequences_padded.npy')

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

# Construct 
RNJ_tic = time.process_time()
taxa = [str(i) for i in range(n_leaves)]
RNJ_tree = relaxed_neighbor_joining(dist_matrix, taxa)
RNJ_toc = time.process_time()

# Save 
pickle.dump(RNJ_tree, open('results/RNJ_tree.pickle', 'wb'))

# Print 
print('*******************')
print('DONE')
print('Distance matrix time= ', way_1_toc - way_1_tic)
print('RNJ clustering time= ', RNJ_toc - RNJ_tic)
print('*******************')
